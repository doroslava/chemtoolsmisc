from e3fp.pipeline import fprints_from_sdf
from os import listdir
from os.path import isfile
from os.path import join
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit import Chem
import pandas
from e3fp.pipeline import fprints_from_smiles
from glob import glob
from python_utilities.parallel import Parallelizer
from e3fp.conformer.util import smiles_to_dict
from python_utilities.parallel import make_data_iterator


#parameters for fingerprint generation
fprint_params = {'bits': 1024, 'rdkit_invariants': True, "first": 25}
confgen_params = {"num_conf": 50}

#Creating fingerprints from the modulators from the crystal structure
modulatorSDFSet = [f for f in listdir("/home/dsribar/Documents/tlr8/Repurposing/Conformers_CS") if
             isfile(join("/home/dsribar/Documents/tlr8/Repurposing/Conformers_CS", f))]

modulatorFpsList= []

for f in modulatorSDFSet:
    fp = fprints_from_sdf("/home/dsribar/Documents/tlr8/Repurposing/Conformers_CS/"+f,fprint_params=fprint_params)
    modulatorFpsList.append(fp)

modulatorFps = []
modulatorIDs = []
for sublist in modulatorFpsList:
    for item in sublist:
        modulatorFps.append(item.to_rdkit())
        modulatorIDs.append(item.name)


#rdkit diversity picker
#nfps=len(modulatorFps)

#def distij(i,j,fps=modulatorFps):
#    return 1-DataStructs.TanimotoSimilarity(fps[i],fps[j])
#picker = MaxMinPicker()
#pickIndices = list(picker.LazyPick(distij,nfps,7,seed=23))

#pickedModulatorIDs=[modulatorIDs[x] for x in pickIndices]

#print "Number of different conformations from the crystal structures : ", nfps

#print "Most diverse conformations are: ",pickedModulatorIDs, "with indices ", pickIndices

#########
# Read the sdf file downloaded from DrugBank (the biotech drugs are already excluded), write the new .smi file with smiles
# and corresponding identifier for each molecule in the database
#########
print "Reading the sdf data set..."

librarySDFSet = [x for x in Chem.SDMolSupplier('/mdspace/dora/TLR8/Old_data/Repurposing/Drug_bank/structures.sdf')]
while librarySDFSet.count(None): librarySDFSet.remove(None)

smiles = [x.GetProp("SMILES") for x in librarySDFSet]
drugbankID = [x.GetProp("DRUGBANK_ID") for x in librarySDFSet]


smiDataFrame = pandas.DataFrame(
        {"col1":smiles,
        "col2": drugbankID})

smiDataFrame.to_csv("/mdspace/dora/TLR8/Old_data/Repurposing/Drug_bank/drugbank_structures.smi",index=False,sep="\t",header=False)

##########
# Generating conformations and calculating 3D Fingerprints for the molecules in the smi.format, generated in the previous step
#There are still some molecules with more substantially rotable bonds, which slows the conformation generation time. Therefore,
# only molecules with smaller number of rotable bond than treshold are considered
##########
rotableBonds = [int(x.GetProp("JCHEM_ROTATABLE_BOND_COUNT")) for x in librarySDFSet]

rotableBondTreshold = 15
sdfSetSize = len(librarySDFSet)

print "Generating conformations and 3D fingerprints..."
print "Treshold for rotable bonds: ",rotableBondTreshold
print "Size of the sdf set: ", sdfSetSize

libraryFpsList = []
FailedFpsList = []




for x in range(0,sdfSetSize):
    if rotableBonds[x]<rotableBondTreshold:
        try:
            libraryFpsList.append(fprints_from_smiles(smiles[x], drugbankID[x], confgen_params=confgen_params, fprint_params=fprint_params, save=False))
        except:
            FailedFpsList.append(x)
    else:
        FailedFpsList.append(x)


smiles = []
drugbankID =[]
drugbankName = []
drugGroup = []


for i,x in enumerate(librarySDFSet):
        if i not in FailedFpsList:
            drugbankID.append(x.GetProp("DRUGBANK_ID"))
            smiles.append(x.GetProp("SMILES"))
            drugbankName.append(x.GetProp("GENERIC_NAME"))
            drugGroup.append(x.GetProp("DRUG_GROUPS"))


libraryFps = []
libraryFpsIDs = []

for sublist in libraryFpsList:
    for item in sublist:
        libraryFps.append(item.to_rdkit())
        libraryFpsIDs.append(item.name)


#Similiarity searching

drugbankDataFrame = pandas.read_csv("/home/dsribar/Documents/tlr8/Repurposing/drug_target.csv")

libraryFpsIDstidy = []

for f in libraryFpsIDs:
    libraryFpsIDstidy.append(f.split("_")[0])

for i in range(0,len(modulatorFps)):
    similiarity=[]
    for x in range(0,len(libraryFps)):
        similiarity.append(DataStructs.TanimotoSimilarity(modulatorFps[i],libraryFps[x]))
    similiarityDataFrame = pandas.DataFrame(
        {'DrugBank ID': libraryFpsIDstidy,'Tanimoto similiarity': similiarity})
#   print similiarityDataFrame
    similiarityDataFrameGrouped=similiarityDataFrame.groupby(['DrugBank ID'], sort=False)["DrugBank ID","Tanimoto similiarity"].max()
    similiarityDataFrameGrouped["Generic Name"] = drugbankName
    similiarityDataFrameGrouped["Drug Group"] = drugGroup
    similiarityDataFrameGrouped["Smiles"] = smiles
    finalDataFrame = pandas.merge(similiarityDataFrameGrouped, drugbankDataFrame, on="DrugBank ID")
    finalDataFrameTopHits = finalDataFrame.sort_values(by="Tanimoto similiarity", ascending=False)
    finalDataFrameTopHits = finalDataFrameTopHits[["DrugBank ID", "Generic Name", "Tanimoto similiarity", "Drug Group", "UniProt Name","Smiles"]]
    finalDataFrameTopHits.to_csv("/home/dsribar/Documents/tlr8/Repurposing/3D_hits/50_confs_25_best/"+modulatorIDs[i]+".csv",index=False)
