import subprocess
import pandas as pd


def glasso(otu_tables,lam = 0.03,normtype = "dirichlet",chunkyness = 0,name = "glasso",table_type="counts",return_precision = False):

    """Python function to run normalCowboy from script.

    :param otu_tables: set of otu tables (if multiple data types) or single table. (read counts or relative abundance)
    :type otu_tables: dict(pd.DataFrame) or pd.DataFrame

    :param lam: LASSO penalty. Default 0.03
    :type lam: float

    :param normtype: type of table normalization. Options - dirichlet/add_pseudocounts/robust. Default "dirichlet"
    :type normtype: str

    :param chunkyness: Size of sub-blocks to fit. 0 for regular GLASSO. Default 0
    :type chunkyness: int

    :param name: name to use for temp files passed to julia. Default "glasso"
    :type name: str

    :param table_type: type of OTU table(s) provided. Options counts/relative. Default "counts"
    :type table_type: str

    :param return_precision:  Whether to return the precision matrix instead of the covariance matrix. Default False
    :type return_precision: bool

    :return: estimated covariance matrix or precision matrix
    :rtype: pd.DataFrame
    """

    tot_taxa = 0

    if return_precision:
        rtype = "precision"
    else:
        rtype = "covariance"

    if isinstance(otu_tables,dict):
        csvnames = []
        for ky,tab in otu_tables.items():
            tot_taxa += len(tab)
            tab.to_csv("otutab_{}.csv".format(ky))
            csvnames += ["otutab_{}.csv".format(ky)]
        command = "julia ../normalcowboy.jl/normalCowboy.jl " + " ".join(csvnames) + " -l={} -c={} -z={} -n={} -t={} -r={}".format(lam,chunkyness,normtype,name,table_type,rtype)
    
    else:
        tot_taxa = len(otu_tables)
        otu_tables.to_csv("otutab.csv")
        command = "julia ../normalcowboy.jl/normalCowboy.jl otutab.csv" + " -l={} -c={} -z={} -n={} -t={} -r={}".format(lam,chunkyness,normtype,name,table_type,rtype)

    print("[glasso] Computing GLASSO Fit using NormalCowboy with {} total taxa".format(tot_taxa))
    subprocess.run(command,shell = True)
    glasso_fit = pd.read_csv("{}.csv".format(name),index_col = 0)
    return glasso_fit

