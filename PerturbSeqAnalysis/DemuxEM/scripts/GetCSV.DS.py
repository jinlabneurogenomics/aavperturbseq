from pegasusio import read_input
import sys
import pandas as pd
from numpy.random import binomial

args=sys.argv
RNAh5=args[1]
outCSV=args[2]
numToDS=int(args[3])

adata=read_input(RNAh5,select_modality="crispr").to_anndata()
counts=adata.X.sum(axis=0).tolist()[0] ##get number of UMI per guide
adata=adata[:,[c>0 for c in counts]] ##exclude ones with no counts
counts=[c for c in counts if c>0] ##get number of UMI per guide for guides with >0 UMI
#numToDS=min(counts)
DS=[min(numToDS/c,1) for c in counts] ##The amound to downsample each guide by
X=adata.X.todense()
print("Downsample")
X=binomial(X,DS)
print("Final counts")
print(X.sum(axis=0).tolist()[0])

print("Save")
dat=pd.DataFrame(X.T,columns=adata.obs.index,index=[i for i in adata.var.index]).reset_index()
dat.to_csv(outCSV,index=False)