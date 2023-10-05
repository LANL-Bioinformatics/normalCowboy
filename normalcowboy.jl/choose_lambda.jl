using ArgParse
using CSV,Tables

include("model_selection.jl")


s = ArgParseSettings()
@add_arg_table s begin
    "--chunksize", "-c"
        help = "Size of sub-blocks to fit"
        default = 0
        arg_type = Int64
    "--zero_method","-z"
        help = "How zeros are handled - options dirichlet/add_pseudocounts/robust"
        default = "dirichlet"
    "--criterion", "-m"
        help = "Model selection criterion (method): StARS, BIC, or AIC"
        default = "StARS"
    "--StARSN", "-N"
        help = "Subsample number for StARS (Only used if criterion is StARS)"
        default = 50
        arg_type = Int64
    "--StARSb", "-b"
        help = "Subsample size for StARS, as percentage of total number of samples (Only used if criterion is StARS)"
        default = 0.8
        arg_type = Float64
    "--Subtable","-s"
        help = "Number of nodes to include, choose 0 for all nodes"
        default = 0
        arg_type = Int64
    "--Repeat","-r"
        help = "Number of times to repeat model selection. Returns average."
        default = 1
        arg_type = Int64
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

    criterion = parsed_args["criterion"]
    zmethod = parsed_args["zero_method"]
    subt = parsed_args["Subtable"]
    SN = parsed_args["StARSN"]
    Sb = parsed_args["StARSb"]
    chsize = parsed_args["chunksize"]

    tlam = 0
    for i = 1:parsed_args["Repeat"]
        tlam = tlam + model_selection(dsets,criterion;method = zmethod,S = subt,N = SN,Stb = Sb,beta = 0.05,r = 1000,minl = 0.01,chunkyness=chsize)
    end
    alam = tlam/parsed_args["Repeat"]
    println(alam)

end


main()