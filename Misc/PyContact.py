from PyContact.core.Scripting import PyContactJob, JobConfig
import MDAnalysis as mda
import os
import pandas
import seaborn as sb
import copy
import matplotlib.pyplot as plt

#Before running the script be sure to both activate environment in command line and import it in pycharm


### Definition of functions for running th pycontact


#Function which includes pycontact scripting, input the folder with replicates, and returns corresponding dataframe with
#bond occurance for replicates


def pycontact(folder,replicates=3):
    dataFinal = pandas.DataFrame()
    dataFrame = {}
    for i in range(0,replicates):
        data = pandas.DataFrame ()
        structureFile=folder+"pycontact"+str(i)+".pdb"
        trajectoryFile=folder+"pycontact"+str(i)+".dcd"
        job = PyContactJob(structureFile,trajectoryFile, "title",
                       JobConfig(5.0, 3.5, 120, [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], "segid A", "segid B"))
        job.runJob(1)
        contacts = job.analyzer.finalAccumulatedContacts
        for x in contacts:
            if (x.contactTypeAsShortcut()!="other"):
                if (x.contactTypeAsShortcut!="hbond"):
                    occurrence = float(len(x.getScoreArray()) - x.getScoreArray().count(0)) / len(x.getScoreArray()) * 100
                else:
                    occurrence=x.hbond_percentage()
                pair= x.human_readable_title() + "--" + str(x.determine_ctype())
                data = data.append(pandas.Series([x.contactTypeAsShortcut(),pair,occurrence]), ignore_index=True)

        dataFrame[i]=data
        job = None

    dataFinal =  dataFrame[0].merge(dataFrame[1].merge(dataFrame[2], on=[0,1], how="outer"),on=[0,1], how="outer")
    dataFinal.columns = ["bondtype","pair","replicate1","replicate2", "replicate3"]
    return(dataFinal)
#Filters the occurence dataframe based on the occurence in replicates
def getCommonPairs(HBondPair, treshold=20):
    return HBondPair[(HBondPair.replicate1>treshold)&(HBondPair.replicate2>treshold)&(HBondPair.replicate3>treshold)]

#Plot the occurence
def plotPairs(df,bondtype):
    newdf=[]
    for i in range(0,len(df)):
        newdf.append(pandas.DataFrame(df[i].loc[df[i]["bondtype"]==bondtype].values.T[2:], columns=df[i].loc[df[i]["bondtype"]==bondtype].pair.tolist()))
        newdf[i]["category"]=i+1
    dataFrame = reduce(lambda left, right: pandas.merge(left, right,how='outer'), newdf)
    tmp= dataFrame.pop("category")
    dataFrame["category"]=tmp
    if len(dataFrame.columns)<20:
        ncol,nrow=10,2
    else:
        ncol,nrow=10,3
    f, axes = plt.subplots(nrow,ncol,sharey=True)
    br = 0
    for i in range(0,len(dataFrame.columns)-1):
        if ((i%30 == 0) and (i/30>0)):
            f, axes = plt.subplots(nrow,ncol, sharey=True)
            br+=30
        sb.stripplot(y=dataFrame.columns[i], data=dataFrame, x="category", jitter=True, ax=axes[(i-br)/10, (i-br) % 10])
    return(dataFrame)




bondsApo= pycontact("/mdspace/dsribar/TLR8/MD/DESMOND/Apo/Apo/Analysis/pycontact/")

bondsCUCPT9= pycontact("/mdspace/dsribar/TLR8/MD/DESMOND/CUCPT9/Analysis/pycontact/")

bondsCL097=pycontact("/mdspace/dsribar/TLR8/MD/DESMOND/CL097/CL097/Analysis/pycontact/")

bondsHolo=pycontact("/mdspace/dsribar/TLR8/MD/DESMOND/CL097/Holo/Analysis/pycontact/")

bondsssRNA=pycontact("/mdspace/dsribar/TLR8/MD/DESMOND/ssRNA/noions/Analysis/pycontact/")



daf=[getCommonPairs(bondsApo),getCommonPairs(bondsCUCPT9),getCommonPairs(bondsCL097),getCommonPairs(bondsHolo),getCommonPairs(bondsssRNA)]

daf=[bondsApo,bondsCUCPT9,bondsCL097,bondsHolo,bondsssRNA]

plotPairs(copy.deepcopy(daf),"hbond")
plt.show()








