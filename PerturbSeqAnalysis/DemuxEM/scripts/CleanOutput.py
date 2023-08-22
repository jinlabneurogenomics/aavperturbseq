from pegasusio import read_input
import sys
import pandas as pd

args=sys.argv
inZarr=args[1]
outCSV=args[2]

dat=read_input(inZarr).obs.reset_index()
dat.to_csv(outCSV,index=False,sep="\t")