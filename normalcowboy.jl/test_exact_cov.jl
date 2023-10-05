using ArgParse
using CSV,Tables

include("makeTable.jl")
include("solveProgram.jl")


#######        How do mean/variance of diagonal chunks change with chunk size

s = ArgParseSettings()
@add_arg_table s begin
    "--lasso", "-l"
        help = "LASSO penalty strength."
        default = 1.0
        arg_type = Float64
    "--chunksize", "-c"
        help = "Size of sub-blocks to fit"
        default = 10
        arg_type = Int64
    "--name","-n"
        help = "Name of output network"
        default = "network"
    "covariance"
        help = "exact covariances to fit to network."
        required = true
end
parsed_args = parse_args(ARGS, s)

function import_covariance(TableFile::String)
    #=
    read in covariance table .csv file to Matrix.
    =#
    csv_reader = CSV.File(TableFile)
    data = DataFrame(csv_reader)
    datamat = Matrix(data[:,2:end])
    rwnms = string.(data[:,1])
    return datamat,rwnms
end

function main()

    cov_table_pth = parsed_args["covariance"]

    cov_table,nms = import_covariance(cov_table_pth)
    
    if parsed_args["chunksize"] !=0
        # precision_m = solve_CholChunky(cov_table,parsed_args["lasso"],chnksize = parsed_args["chunksize"])
        precision_m = solve_Chunky(cov_table,parsed_args["lasso"],chnksize = parsed_args["chunksize"])


    else
        # precision_m = solve_CholClean(cov_table,parsed_args["lasso"])# solve_withSDP(cov_table,parsed_args["lasso"])#
        precision_m =  solve_withSDP(cov_table,parsed_args["lasso"])
    end

    nm = parsed_args["name"]

    cov_table_tab = Tables.table(hcat(nms,cov_table))
    CSV.write("$nm.cov_table.csv", cov_table_tab, header = vcat([""],nms))

    precision_table = Tables.table(hcat(nms,precision_m))
    CSV.write("$nm.csv",  precision_table, header = vcat([""],nms))

    # marg_table = Tables.table(hcat(clr_table[:otus],off_diag))
    # CSV.write("$nm.marginals.csv",marg_table,header = vcat([""],clr_table[:otus]))

end


main()
