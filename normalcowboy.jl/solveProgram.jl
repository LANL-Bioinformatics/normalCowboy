using JuMP
using SCS
import MathOptInterface as MOI
# using Gurobi
using LinearAlgebra
using BenchmarkTools
using ProgressMeter
using Combinatorics
using Alpine
using EAGO


function solve_withSDP(M::Matrix,l::Number;solver = SCS.Optimizer,bequiet = true)

    model = Model(solver)
    if bequiet
        set_silent(model)
    end
    sz = size(M)[1]

    # println(size(M))

    @variable(model, X[1:sz,1:sz],PSD)
    @constraint(model,[0;1;vec(X)] in MOI.LogDetConeSquare(sz))#constrain log-det > 1
    @variable(model, t)
    @constraint(model, [t; vec(X)] in MOI.NormOneCone(1 + length(vec(X))))#t >= sum(abs(X))
    @objective(model,Min, tr(X*M) + l*t)

    # println("OPTIMIZING!")

    optimize!(model)

    Sol = value(X)

    return (sz/tr(Sol*M))*Sol

end

function innerloop_SDP!(ij,M,l,chunks,marginal,diag_blocks,diag_counter,solver,bequiet)
    i = ij[1]
    j = ij[2]
    indx = union(chunks[i],chunks[j])
    blk = solve_withSDP(M[indx,indx],l,solver=solver,bequiet = bequiet)
    # all_fits[union(chunks[i],chunks[j])] = blk

    marginal[indx,indx] = inv(blk)

    diag_blocks[chunks[i]][:,:,diag_counter[chunks[i]]] = marginal[chunks[i],chunks[i]]
    diag_blocks[chunks[j]][:,:,diag_counter[chunks[j]]] = marginal[chunks[j],chunks[j]]
    diag_counter[chunks[i]] += 1
    diag_counter[chunks[j]] += 1
end

function solve_Chunky(M,l;chnksize = 10,solver=SCS.Optimizer,bequiet = true)

    numrngs = Int(floor(size(M)[1]/chnksize))

    chunks = [(i-1)*chnksize+1:i*chnksize for i in 1:numrngs]
    nmCh = length(chunks)
    if chnksize*numrngs < size(M)[1]
        append!(chunks,[range(chnksize*numrngs + 1,size(M)[1])])
    end

    marginal = zeros(size(M))
    diag_blocks = Dict(chunks[i]=>zeros(chnksize,chnksize,nmCh-1) for i in 1:nmCh)
    diag_counter = Dict(chunks[i]=>1 for i in 1:nmCh)


    @showprogress for ij in combinations(1:nmCh,2)
        innerloop_SDP!(ij,M,l,chunks,marginal,diag_blocks,diag_counter,solver,bequiet)
    end

    the_net = copy(marginal)#/(nmCh-1)
    for ky in keys(diag_blocks)
        the_net[ky,ky] = mean(diag_blocks[ky],dims = 3)[:,:,1]
    end

    return the_net#,diag_blocks,marginal

end

# function solve_withQuadratic(M::Matrix,l::Number;solver=Gurobi.Optimizer)

#     model = Model(solver)

#     sz = size(M)[1]

#     @variable(model,X[1:sz,1:sz],Symmetric)
#     @constraint(model,diagc[i=1:sz], X[i,i] >= 0.1 )
#     @constraint(model,quadc[i=1:sz,j=1:sz], X[i,i]*X[j,j] >= X[i,j]^2)
#     ### Do I still need the full rank constraint?

#     @constraint(model,[0;1;vec(X)] in MOI.LogDetConeSquare(sz))#Are there more efficient ways to constrain rank? Probably!
    

#     @objective(model,Min,tr(X*M) + l*sum(X))

#     optimize!(model)

#     return value(X)

# end

function solve_withCholesky(M::Matrix,l::Number;solver=EAGO.Optimizer,bequiet = true)

    model = Model(solver)

    if bequiet
        set_silent(model)
    end

    sz = size(M)[1]

    @variable(model,A[1:sz,1:sz],Symmetric)
    C = ones(sz,sz).*LowerTriangular(A)

    @constraint(model,diagc[i=1:sz],C[i,i]>=0)
    # @variable(model,diags[1:sz])
    # @NLconstraint(model,logeq[i=1:sz],diags[i] == log(C[i,i]))
    # @constraint(model,det,sum(diags[j] for j=1:sz) >= 1)

    @NLconstraint(model,det,prod(C[j,j] for j =1:sz) >= 1)
    @objective(model,Min,tr((C*transpose(C))*M)+l*sum(C*transpose(C)))

    optimize!(model)

    CV = value.(C)

    X = CV*transpose(CV)

    # println(tr(X*M))
    # println(det(X))


    for i = 1:sz
        println((sqrt(sz/tr(X*M)))*CV[i,:])
    end

    return (sz/tr(X*M))*X

end

# function solve_QuadClean(M,l)

#     X = solve_withQuadratic(M,l,solver = () -> Gurobi.Optimizer(env))

#     return X

# end

# function solve_CholClean(M,l,env)

#     X = solve_withCholesky(M,l,solver = () -> Gurobi.Optimizer(env))

#     return X

# end

function innerloop_Chol!(ij,M,l,chunks,marginal,diag_blocks,diag_counter,solver,bequiet)
    i = ij[1]
    j = ij[2]
    blk = solve_withCholesky(M[union(chunks[i],chunks[j]),union(chunks[i],chunks[j])],l,solver=solver,bequiet = bequiet)

    marginal[union(chunks[i],chunks[j]),union(chunks[i],chunks[j])] = inv(blk)

    diag_blocks[chunks[i]][diag_counter[chunks[i]],:,:] = marginal[chunks[i],chunks[i]]
    diag_blocks[chunks[j]][diag_counter[chunks[j]],:,:] = marginal[chunks[j],chunks[j]]
    diag_counter[chunks[i]] += 1
    diag_counter[chunks[j]] += 1
end


function solve_CholChunky(M,l,env;chnksize = 10,solver = EAGO.Optimizer,bequiet=true)

    numrngs = Int(floor(size(M)[1]/chnksize))

    chunks = [(i-1)*chnksize+1:i*chnksize for i in 1:numrngs]
    nmCh = length(chunks)
    if chnksize*numrngs < size(M)[1]
        append!(chunks,[range(chnksize*numrngs + 1,size(M)[1])])
    end

    marginal = zeros(size(M))
    diag_blocks = Dict(chunks[i]=>zeros(nmCh-1,chnksize,chnksize) for i in 1:nmCh)
    diag_counter = Dict(chunks[i]=>1 for i in 1:nmCh)

    @showprogress for ij in combinations(1:nmCh,2)

        innerloop_Chol!(ij,M,l,chunks,marginal,diag_blocks,diag_counter,solver,bequiet)
        # i = ij[1]
        # j = ij[2]
        # blk = solve_CholClean(M[union(chunks[i],chunks[j]),union(chunks[i],chunks[j])],l,env)

        # marginal[union(chunks[i],chunks[j]),union(chunks[i],chunks[j])] = inv(blk)
        # diag_blocks[chunks[i]][diag_counter[chunks[i]],:,:] = marginal[chunks[i],chunks[i]]
        # diag_blocks[chunks[j]][diag_counter[chunks[j]],:,:] = marginal[chunks[j],chunks[j]]
        # diag_counter[chunks[i]] += 1
        # diag_counter[chunks[j]] += 1
    end

    the_net = copy(marginal)#/(nmCh-1)
    for ky in keys(diag_blocks)
        the_net[ky,ky] = mean(diag_blocks[ky],dims = 1)[1,:,:]
    end

    return the_net

end

