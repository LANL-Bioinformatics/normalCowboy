using LinearAlgebra
using Optim
using Roots

include("makeTable.jl")
include("solveProgram.jl")

function compute_AIC(est,emp)
    L = log(det(est)) + tr(inv(est)*emp)
    m = sum(abs.(est).>10^-8)
    p = size(emp)[1]
    d = m*(m-1)/2 + p
    return 2*L + 2*d
end

function est_and_AIC(l,cov_table;chunkyness = 0)
    if chunkyness !=0
        est = solve_Chunky(cov_table,l,chnksize = chunkyness)
    else
        precision_m = solve_withSDP(cov_table,l)
        est = inv(precision_m)
    end
    aic = compute_AIC(est,cov_table)
    return aic
end

function choose_AIC(cov_table;chunkyness = 0)
    res = optimize(l->est_and_AIC(l,cov_table,chunkyness=chunkyness),0,1)
    return Optim.minimizer(res)
end



function compute_BIC(est,emp,n)
    L = log(det(est)) + tr(inv(est)*emp)
    m = sum(abs.(est).>10^-8)
    p = size(emp)[1]
    d = m*(m-1)/2 + p
    return 2*L + log(n)*d
end

function est_and_BIC(l,cov_table,n;chunkyness = 0)
    if chunkyness !=0
        est = solve_Chunky(cov_table,l,chnksize = chunkyness)
    else
        precision_m = solve_withSDP(cov_table,l)
        est = inv(precision_m)
    end
    bic = compute_BIC(est,cov_table,n)
    return bic
end

function choose_BIC(cov_table,n;chunkyness = 0)
    res = optimize(l->est_and_BIC(l,cov_table,n,chunkyness=chunkyness),0,1)
    return Optim.minimizer(res)
end



function do_fit(clr_tab,l;chunkyness = 0)
    cov_tab = cov(clr_tab,dims=2)
    if chunkyness !=0
        est = solve_Chunky(cov_tab,l,chnksize = chunkyness)
    else
        precision_m = solve_withSDP(cov_tab,l)
        est = inv(precision_m)
    end
    return est
end

function compute_Dhat(clr_tab,l,N,b;chunkyness = 0,thr = 10^-8)
    p,ns = size(clr_tab)
    nm = zeros(p,p)
    for i in 1:N
        tb = do_fit(clr_tab[:,sample(1:ns,b,replace=false)],l,chunkyness=chunkyness)
        pr = abs.(tb).>thr
        nm = nm + pr
    end
    th = nm ./ N
    xis = 2 .* th .* (1 .- th)
    dhat = sum(xis)/(p*(p-1))
    return dhat
end

function choose_Stars(clr_tab,N,b;beta = 0.05,r = 100,maxl = 0.5,minl = 0.01,chunkyness=0)
    strtl = (maxl+minl)/2
    l = 0
    Dh0 = compute_Dhat(clr_tab,strtl,N,b,chunkyness=chunkyness)
    # println(Dh0)
    # println("======")
    if Dh0 > beta #assume we started lower than the correct l
        # println(strtl,":",r,":",maxl)
        # println(collect(strtl:r:maxl))
        for outer l = strtl:1/r:maxl
            Dh = compute_Dhat(clr_tab,l,N,b,chunkyness=chunkyness)
            if Dh < beta
                break
            end
        end
    end
    if Dh0 < beta #assume we started higher (this could be wrong if maxl is too low!)
        for outer l = reverse(minl:1/r:strtl)
            Dh = compute_Dhat(clr_tab,l,N,b,chunkyness=chunkyness)
            if Dh > beta
                break
            end
        end
    end
    return l
end

# function choose_Stars_Root(clr_tab,N,b;beta = 0.05,r = 100,maxl = 0.5,minl = 0.01,chunkyness=0)
 
# end

function model_selection(dsets,criterion;method = "dirichlet",S = 0,N = 100,Stb = 0.8,beta = 0.05,r = 1000,minl = 0.01,chunkyness=0,maxl = 0.2)
    cov_table,clr_table_dat = make_covtable(dsets,method=method)
    clr_table = clr_table_dat[:otuTable]

    if criterion == "StARS"

        ns = size(clr_table)[2]
        b = Int(floor(Stb*ns))

        if S>0
            nt = size(clr_table)[1]
            clr_table = clr_table[sample(1:nt,S,replace=false),:]
        end

        lam = choose_Stars(clr_table,N,b;beta = beta, r = r,minl = minl,maxl=maxl,chunkyness=chunkyness)

        return lam
    
    elseif criterion == "AIC"

        if S>0
            nt = size(cov_table)[1]
            smpl = sample(1:nt,S,replace=false)
            cov_table = cov_table[smpl,smpl]
        end

        lam = choose_AIC(cov_table,chunkyness=chunkyness)
        return lam
        
    elseif criterion == "BIC"

        if S>0
            nt = size(cov_table)[1]
            smpl = sample(1:nt,S,replace=false)
            cov_table = cov_table[smpl,smpl]
        end
        
        ns = size(clr_table)[2]
        lam = choose_BIC(cov_table,ns,chunkyness=chunkyness)
        return lam

    end
end