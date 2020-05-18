from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
import pandas
from os import listdir
from os.path import isfile
from os.path import join

fingerprintType="3D"

onlyFiles = [f for f in listdir("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/50_confs_25_best") if
             isfile(join("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/50_confs_25_best", f))]

frames = []
queryFile = []
for i in onlyFiles:
    fileToRead = join("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/50_confs_25_best",i)
    DataFrame = pandas.read_csv(fileToRead,nrows=50)
    subset = DataFrame[["DrugBank ID","Tanimoto similiarity","Generic Name","UniProt Name","Smiles"]]
    frames.append(DataFrame)
    queryFile.append([i]*len(DataFrame))


flattenedQueryFile = [y for x in queryFile for y in x]

result = pandas.concat(frames)

result["Query"] = flattenedQueryFile

result["Fingerprint Type"] = [fingerprintType]*len(result)

#resultNodupsTable = result.drop_duplicates(subset="DrugBank ID")

resultNodupsTableSorted= result[result.groupby(['DrugBank ID'])['Tanimoto similiarity'].transform(max) == result['Tanimoto similiarity']]


resultNodupsTableSorted.to_csv("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/50_confs_25_best/50_confs_25_best.csv",index=False,sep="\t",header=True)


DataFrame2D = pandas.read_csv("/home/dsribar/Documents/tlr8/Repurposing/2D_hits/2D_best.csv",sep="\t")
DataFrame50 = pandas.read_csv("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/50_confs_25_best/50_confs_25_best.csv",sep="\t")
DataFrame25 = pandas.read_csv("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/25_confs_5_best/25_confs_5_best.csv",sep="\t")

list(set(DataFrame25).difference(DataFrame5))

a = DataFrame2D["DrugBank ID"].tolist()
b = DataFrame50["DrugBank ID"].tolist()
c = DataFrame25["DrugBank ID"].tolist()

print len(set(b)&set(c))

print len((set(b)& set(c)))+len((set(b).difference(c)))+len((set(c).difference(b)))

unique50=set(b).difference(c)
unique25=set(c).difference(b)

commonHits = list(set(a)& set(c))
UniqueHits2D = list(set(a).difference(c))
UniqueHits3D = list(set(c).difference(a))



groups = DataFrame50["Drug Group"].tolist()
ids =  DataFrame50["DrugBank ID"].tolist()
names = DataFrame50["Generic Name"].tolist()

dataFrame = pandas.DataFrame({"Groups": groups,"id":ids,"Name":names})
dataFrame = dataFrame.loc[(dataFrame.Groups !="experimental")]
dataFrame[dataFrame["id"].isin(unique50)]






commonHitsSmiles = DataFrame2D[DataFrame2D["DrugBank ID"].isin(commonHits)][["Smiles","DrugBank ID"]].drop_duplicates(subset="DrugBank ID")
UniqueHits2DSmiles = DataFrame2D[DataFrame2D["DrugBank ID"].isin(UniqueHits2D)][["Smiles","DrugBank ID"]].drop_duplicates(subset="DrugBank ID")
UniqueHits3DSmiles = DataFrame25[DataFrame25["DrugBank ID"].isin(UniqueHits3D)][["Smiles","DrugBank ID"]].drop_duplicates(subset="DrugBank ID")

commonHitsSmiles.to_csv("/home/dsribar/Documents/tlr8/Repurposing/commonHits.smi",index=False,sep="\t",header=False)
UniqueHits2DSmiles.to_csv("/home/dsribar/Documents/tlr8/Repurposing/UniqueHits2D.smi",index=False,sep="\t",header=False)
UniqueHits3DSmiles.to_csv("/home/dsribar/Documents/tlr8/Repurposing/UniqueHits3D.smi",index=False,sep="\t",header=False)
