import pysam
import pandas as pd
import os

##
##Gets insertion and deletions for all bams in directory, excluding those that aren't from perturbations
##Use "../Ch[1-5]" as dir
##
def GetInsertDel_directory(dir,chrom,start,end):
    os.system("ls -1 "+dir+"/*bam | grep -v reads | grep -v combined > temp.txt")
    bams=open("temp.txt","r").read().strip().split("\n") ##Get list of bams
    print(bams)
    results=[GetInsertDel(bamfile,chrom,start,end) for bamfile in bams]
    for i in range(0,len(bams)):
        results[i]["Name"]=["Ch"+str(i+1)+"_"+nam for nam in results[i]["CBC"]]
    dat=pd.concat(results,axis=0)
    return(dat)





##
##Runs GetInsertDel_directory for multiple regions
##
def GetInsertDel_regions(dir,regs):
    results=[GetInsertDel_directory(dir,reg[0],reg[1],reg[2]) for reg in regs]
    for i in range(0,len(regs)):
        reg=regs[i]
        results[i]["Region"]=reg[0]+":"+str(reg[1])+"-"+str(reg[2])
        results[i]["Guide"]=reg[3]
    dat=pd.concat(results,axis=0)
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
    ret=[]
    for seq in samfile.fetch(contig=chrom,start=start,end=end):
        if not seq.has_tag("CB"):
            continue;
        cbc=seq.get_tag("CB")
        [insert,delete,overlaps]=GetInsertDelete(seq,chrom,start,end)
        if overlaps==0:
            continue;
        numTot=numTot+1
        numInsert=numInsert+insert
        numDelete=numDelete+delete
        if insert>0 and delete>0:
            numBoth=numBoth+1
        ret.append([insert,delete,cbc])
    samfile.close()
    print(bamfile)
    print("Total "+str(numTot))
    print("Insert "+str(numInsert))
    print("Delete "+str(numDelete))
    print(" ")
    dat=pd.DataFrame(ret)
    dat.columns=["Insert","Delete","CBC"]
    dat["Reads"]=1
    dat=dat.groupby("CBC").sum().reset_index()
    return(dat)

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

        
                
        


