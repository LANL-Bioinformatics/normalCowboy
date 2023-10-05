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
        help = "Size of sub-blocks to fit. 0 for regular GLASSO"
        default = 0
        arg_type = Int64
    "--name","-n"
        help = "Name of output network"
        default = "network"
    "--zero_method","-z"
        help = "How zeros are handled - options dirichlet/add_pseudocounts/robust"
        default = "dirichlet"
    "--table_type","-t"
        help = "Type of table given - options counts/relative"
        default = "counts"
    "--return_type","-r"
        help = "Type of table returned - covariance or precision (inverse covariance) - options covariance/precision"
        default = "covariance"
    "otutables"
        help = "list of OTU tables used to create network."
        required = true
        nargs='+'
end
parsed_args = parse_args(ARGS, s)

# const env = Gurobi.Env()


function main()
    dsets = Dict()
    for j in parsed_args["otutables"]
        dsets[j] = import_table(j)
    end

    cov_table,clr_table = make_covtable(dsets,method=parsed_args["zero_method"],table_type=parsed_args["table_type"])
    
    if parsed_args["chunksize"] !=0
        covar_matrix = solve_Chunky(cov_table,parsed_args["lasso"],chnksize = parsed_args["chunksize"])
        if lowercase(parsed_args["return_type"]) == "precision"
            precision_m = inv(covar_matrix)
        end
    else
        precision_m = solve_withSDP(cov_table,parsed_args["lasso"])
        if lowercase(parsed_args["return_type"]) == "covariance"
            covar_matrix = inv(precision_m)
        end
    end

    nm = parsed_args["name"]


    cov_table_tab = Tables.table(hcat(clr_table[:otus],cov_table))
    CSV.write("$nm.cov_table.csv", cov_table_tab, header = vcat([""],clr_table[:otus]))

    clr_tab = Tables.table(hcat(clr_table[:otus],clr_table[:otuTable]))
    CSV.write("$nm.clr_table.csv",  clr_tab, header = vcat([""],clr_table[:samples]))

    if lowercase(parsed_args["return_type"]) == "covariance"
        estimate_covtable = Tables.table(hcat(clr_table[:otus],covar_matrix))
        CSV.write("$nm.csv",  estimate_covtable, header = vcat([""],clr_table[:otus]))
    else
        estimate_precision = Tables.table(hcat(clr_table[:otus],precision_m))
        CSV.write("$nm.csv",  estimate_precision, header = vcat([""],clr_table[:otus]))
    end

    # marg_table = Tables.table(hcat(clr_table[:otus],off_diag))
    # CSV.write("$nm.marginals.csv",marg_table,header = vcat([""],clr_table[:otus]))

end


main()
