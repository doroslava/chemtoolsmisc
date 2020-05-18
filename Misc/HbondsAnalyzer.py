#Analyzes data from Hbond analyzer from vmd
def reFormatPycontactInput(path,rep=3,begin1=1,end1=1,begin2=1,end2=1,name1="A",name2="B"):
    os.mkdir(path + "pycontact/")
    for i in range (0,rep):
        structureFile=path+"test"+str(i)+".psf"
        trajectoryFile=path+"test"+str(i)+"_out.trr"

        u=mda.Universe(structureFile,trajectoryFile)

        u.atoms[begin1:end1].segids=name1
        u.atoms[begin2:end2].segids=name2

        A = u.selectAtoms("not resname SPC")

        A.write(path+"pycontact/pycontact"+str(i)+".pdb")

        with mda.Writer(path+"pycontact/pycontact"+str(i)+".dcd", A.n_atoms) as W:
            for ts in u.trajectory:
                W.write(A)
import pandas
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import copy


#Read the data for replicates and  save in one dictionary. Assumes that two files exist-one for number of hydrogen bonds counts for each frame
#in trajectory and one for ocuppancy of each hydrogen bond pair
def stringPercentToFloat(stringPercent):
    return(float(stringPercent.strip('%')) /100)

def readHBondTraj(pathtoFolder):
    HBondTraj={}
    for i in range(3):
        name = 'hbonds' + str(i) + ".dat"
        HBondTraj[i]=pandas.read_csv(os.path.join(pathtoFolder,name), header=None, sep='\s+')
    return(HBondTraj)

def readHBondPairs(pathtoFolder):
    HBondPairs={}
    for i in range(3):
        name = 'hbond-details' + str(i) + ".dat"
        HBondPairs[i]=pandas.read_csv(os.path.join(pathtoFolder,name),sep='\s+',skiprows=1,header=0)
        HBondPairs[i]["occupancy"]=HBondPairs[i]["occupancy"].map(stringPercentToFloat)
        HBondPairs[i]=HBondPairs[i].sort_values(by=["occupancy"],ascending=False)
    dataFrame = HBondPairs[0].merge(HBondPairs[1].merge(HBondPairs[2], on=["donor", "acceptor"], how="outer"),on=["donor", "acceptor"], how="outer")
    dataFrame.columns = ["donor","acceptor","replicate1","replicate2","replicate3"]
    dataFrame.insert(0,"pairs",dataFrame["donor"]+" "+ dataFrame["acceptor"])
    dataFrame= dataFrame.drop(["donor","acceptor"],axis=1)
    return dataFrame

def getCommonPairs(HBondPair, treshold=0.01):
    return HBondPair[(HBondPair.replicate1>treshold)|(HBondPair.replicate2>treshold)| (HBondPair.replicate3>treshold)]


def plotPairs(df):
    for i in range(0,len(df)):
        df[i]=pandas.DataFrame(df[i].values.T[1:], columns=df[i].pairs.tolist())
        df[i]["category"]=i+1
    dataFrame = reduce(lambda left, right: pandas.merge(left, right,how='outer'), df)
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

HBondTrajApo=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/Apo/Analysis")
HBondPairsApo=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/Apo/Analysis")

HBondTrajCUCPT9=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/CUCPT9/Analysis")
HBondPairsCUCPT9=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/CUCPT9/Analysis")

HBondTrajApoCUCPT9=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/CUCPT9/Analysis")
HBondPairsApoCUCPT9=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/CUCPT9/Analysis")

HBondTrajCL0971=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/CL097/pose_1/Analysis")
HBondPairsCL0971=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/CL097/pose_1/Analysis")

HBondTrajCL0972=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/CL097/pose_2/Analysis")
HBondPairsCL0972=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/CL097/pose_2/Analysis")

HBondTrajCL0973=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/CL097/pose_3/Analysis")
HBondPairsCL0973=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/CL097/pose_3/Analysis")

HBondTrajPyrimidine1=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/Pyrimdineagonist/Pose_1/Analysis")
HBondPairsPyrimidine1=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/Pyrimdineagonist/Pose_1/Analysis")

HBondTrajPyrimidine2=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/Apo/Pyrimdineagonist/Pose_2/Analysis")
HBondPairsPyrimidine2=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/Pyrimdineagonist/Pose_2/Analysis")

HBondTrajPyrimidine3=readHBondTraj("//mdspace/dora/TLR8/MD/DESMOND/Apo/Pyrimdineagonist/Pose_3/Analysis")
HBondPairsPyrimidine3=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/Apo/Pyrimdineagonist/Pose_3/Analysis")

HBondTrajCL097=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/CL097/CL097/Analysis")
HBondPairsCL097=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/CL097/CL097/Analysis")

HBondTrajHolo=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/CL097/Holo/Analysis")
HBondPairsHolo=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/CL097/Holo/Analysis")


HBondTrajssRNA=readHBondTraj("/mdspace/dora/TLR8/MD/DESMOND/ssRNA/noions/Analysis")
HBondPairsssRNA=readHBondPairs("/mdspace/dora/TLR8/MD/DESMOND/ssRNA/noions/Analysis")

daf=[data]


plotPairs(copy.deepcopy(daf))



#Example for plotting average number of hydrogen bonds through MD simulation
fig, ax = plt.subplots()

ax.plot(HBondTrajApo[0][0],HBondTrajApo[0][1],".",color="blue")
ax.plot(HBondTrajApo[1][0],HBondTrajApo[1][1],".",color="blue")
ax.plot(HBondTrajApo[2][0],HBondTrajApo[2][1],".",color="blue")

ax.plot(HBondTrajCUCPT9[0][0],HBondTrajCUCPT9[0][1],".",color="yellow")
ax.plot(HBondTrajCUCPT9[1][0],HBondTrajCUCPT9[1][1],".",color="yellow")
ax.plot(HBondTrajCUCPT9[2][0],HBondTrajCUCPT9[2][1],".",color="yellow")







