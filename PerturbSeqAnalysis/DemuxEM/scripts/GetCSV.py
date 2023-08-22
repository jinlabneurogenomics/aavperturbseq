from pegasusio import read_input
import sys
import pandas as pd

args=sys.argv
RNAh5=args[1]
outCSV=args[2]

adata=read_input(RNAh5,select_modality="crispr").to_anndata()
dat=pd.DataFrame(adata.X.todense().T,columns=adata.obs.index,index=[i for i in adata.var.index]).reset_index()
dat.to_csv(outCSV,index=False)