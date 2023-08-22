import pysam
import pandas as pd
import os

##
##Gets insertion and deletions for all bams in directory, excluding those that aren't from perturbations
##
def GetInsertDel_directory(dir,chrom,start,end):
    os.system("ls -1 "+dir+"/*bam | grep -v reads | grep -v combined > temp.txt")
    bams=open("temp.txt","r").read().strip().split("\n") ##Get list of bams
    results=[GetInsertDel(bamfile,chrom,start,end) for bamfile in bams]
    dat=pd.DataFrame(results)
    dat["Bam"]=bams
    dat.columns=["Total","Insert","Delete","Both","Bam"]
    dat["Percent"]=100*((dat["Insert"]+dat["Delete"]-dat["Both"])/dat["Total"])
    dat["dir"]=dir
    dat["chrom"]=chrom
    dat["start"]=start
    dat["end"]=end
    return(dat)

##
##Goes read by read and determines if has insertion/deletion in region of interest
##
def GetInsertDel(bamfile,chrom,start,end):
    samfile = pysam.AlignmentFile(bamfile, "rb") ##Load sam file
    numTot=0 ##Number reads
    numInsert=0 ##Number with insertion
    numDelete=0 ##Number with deletion
    numBoth=0;
    for seq in samfile:
        [insert,delete,overlaps]=GetInsertDelete(seq,chrom,start,end)
        if overlaps==0:
            continue;
        numTot=numTot+1
        numInsert=numInsert+insert
        numDelete=numDelete+delete
        if insert>0 and delete>0:
            numBoth=numBoth+1
    samfile.close()
    print(bamfile)
    print("Total "+str(numTot))
    print("Insert "+str(numInsert))
    print("Delete "+str(numDelete))
    print(" ")
    return([numTot,numInsert,numDelete,numBoth])

##
##Counts number of deletions/insertions based on refPos
##Returns insert (0 or 1), delete (0 or 1), and overlap (0 or 1) where 1 is true, 0 false
##
def GetInsertDelete(seq,chrom,start,end):
    if chrom!=seq.reference_name:
        print("Yuck!")
        print(seq.reference_name)
        return([0,0,0])
    refPos=seq.get_reference_positions(full_length=True) ##the position each base maps to
    insert=0;
    delete=0;
    overlaps=0;
    posprev=None
    pos=None
    inIns=False
    for cur in refPos:
        posprev=pos
        pos=cur
        if pos!=None and pos>start and pos<end:
            if inIns and overlaps==1:
                insert=1;
                inIns=False
            overlaps=1;
            if posprev!=None:
                if pos>posprev+1:
                    delete=1
        if pos==None:
            inIns=True;
        if pos!=None:
            inIns=False
    return([insert,delete,overlaps])

        
                
        


