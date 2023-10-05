using DataFrames
using CSV
using Distributions

function import_table(otuTableFile::String;rows = "OTU")
    #=
    read in otu table .csv file to list of sample names, list of OTUs,
    and Matrix OTU table. Ordering matches between two lists and Matrix.
    Output table is in OTUxSample. Specify input orientation with rows argument.
    =#
    csv_reader = CSV.File(otuTableFile)
    data = DataFrame(csv_reader)
    if lowercase(rows) == "otu"
        otuNames= string.(data[:,1])
        otuTable = data[:,2:end]
        sample_names = names(otuTable)
    elseif lowercase(rows) == "sample"
        sample_names= data[:,1]
        otuTable = data[:,2:end]
        otuNames = names(otuTable)
        otuTable = permutedims(otuTable)
    end
    return (samples=sample_names,otus = otuNames, otuTable=Matrix(otuTable))
end

function split_table(data::NamedTuple,split::Tuple)
    #=
    Splits the data into two seperate data tables, each with the same samples and different OTUs, based on the given split. Split should provide two lists of OTU names, one for each split.
    =#

    #get indices
    set1_ind = findall(otu -> otu in split[1],data[:otus])
    set2_ind = findall(otu -> otu in split[2],data[:otus])

    tab1 = data[:otuTable][set1_ind,:]
    tab2 = data[:otuTable][set2_ind,:]

    # if the ordering of the OTUs as given in split is not the same as in the table, the table will preserve its own ordering. We need to recover that Ordering
    otu1 = data[:otus][set1_ind]
    otu2 = data[:otus][set2_ind]

    return (samples = data[:samples],otus=otu1,otuTable=tab1),(samples=data[:samples],otus=otu2,otuTable=tab2)

end

function combine_table(data1::NamedTuple,data2::NamedTuple)
    #=
    Combines data into a single otu table.
    =#

    allotus = vcat(data1[:otus],data2[:otus])

    if issetequal(data1[:samples],data2[:samples])

        #order of d2 samples in d1
        d2sorder = [findfirst(s -> s == x,data2[:samples]) for x in data1[:samples]]
        combotab = vcat(data1[:otuTable],data2[:otuTable][:,d2sorder])
        return (samples=data1[:samples],otus = allotus,otuTable=combotab)

    else

        println("Removing unmatched samples.")
        allsamples = intersect(data1[:samples],data2[:samples])
        d1sorder = [findfirst(s -> s == x,data1[:samples]) for x in allsamples]
        d2sorder = [findfirst(s -> s == x,data2[:samples]) for x in allsamples]
        combotab = vcat(data1[:otuTable][:,d1sorder],data2[:otuTable][:,d2sorder])
        return (samples=allsamples,otus = allotus,otuTable=combotab)
    end

end

function normalize_table(otuTable::Matrix;method="Dirichlet",table_type="Counts",pseudocount = 1,totcount = 10^5)
    #=
    normalizes otu table and adds psuedocounts. Either adds psuedocount ands normalizes or draws from Dirichlet Distribution using original counts+psuedocount as parameters
    =#

    if lowercase(table_type) == "counts"
        otuTable = otuTable .+ pseudocount
        if lowercase(method) == "dirichlet"
            newtable = float.(copy(otuTable))
            for i in 1:size(newtable)[2]
                newtable[:,i] = rand(Dirichlet(otuTable[:,i]))
            end
            return newtable
        elseif lowercase(method) == "add_pseudocounts"
            newtable = otuTable .+ pseudocount
            return newtable
        else
            return otuTable./sum(otuTable,dims = 1)
        end
    elseif lowercase(table_type) == "relative"
        if lowercase(method) == "dirichlet"
            otuTableCnt = otuTable*(totcount) .+ pseudocount
            newtable = copy(otuTable)
            for i in 1:size(newtable)[2]
                newtable[:,i] = rand(Dirichlet(otuTableCnt[:,i]))
            end
            return newtable
        else
            return otuTable .+ (pseudocount/(totcount))
        end
    end

end

function clr(x::Vector)
    #=
    centered-log-ratio of the vector x
    =#

    log_geo_mean = (1/length(x))*(sum(log.(x)))
    return log.(x) .- log_geo_mean

end

function  robust_clr(x::Vector)
    #=
    "Robust" CLR of the vector x (e.g. ignore 0s)
    =#

    r_log_geo_mean = (1/sum(x.!=0))*(sum(log.(x[x.!=0])))
    rclr = zeros(size(x))
    for i in 1:size(x)[1]
        if x[i] != 0
            rclr[i] = log(x[i])-r_log_geo_mean
        end
    end
    return rclr
    
end

function table_clr(otuTable::Matrix;robust = false)
    #=
    Create table with columns equal to centered-log-ratio of columns of original table.
    =#

    newtable = zeros(size(otuTable))
    for i = 1:size(otuTable)[2]
        if robust
            newtable[:,i] = robust_clr(otuTable[:,i])
        else
            newtable[:,i] = clr(otuTable[:,i])
        end
    end
    return newtable

end

function make_covtable(DataSets::Dict;method="Dirichlet",table_type="Counts",pseudocount = 1,totcount = 10^5)
    #=
    Given a dict of datasets, seperately CLR transforms, reconciles samples, and creates combined covariance table.
    =#

    kylist = collect(keys(DataSets))
    clr_datasets = Dict()

    for ky in kylist
        if lowercase(method) == "robust"
            clr_tab = table_clr(DataSets[ky][:otuTable],robust = true)
        else
            clr_tab = table_clr(normalize_table(DataSets[ky][:otuTable],method=method,table_type=table_type,pseudocount = pseudocount,totcount = totcount))
        end
        clr_datasets[ky] =  (samples=DataSets[ky][:samples],otus = DataSets[ky][:otus], otuTable=clr_tab)
    end

    

    clr_combined = clr_datasets[kylist[1]]
    for i in 2:lastindex(kylist)
        ky = kylist[i]
        clr_combined = combine_table(clr_combined,clr_datasets[ky])
    end

    covtable = cov(clr_combined[:otuTable],dims = 2)

    return covtable,clr_combined

end