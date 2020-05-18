from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
import pandas


#Creating molecules and fingerprints for the library
librarySet1 = [x for x in Chem.SDMolSupplier('/home/dsribar/DIO3_final_selection.sdf')]
while librarySet1.count(None): librarySet1.remove(None)
fps1 = [GetMorganFingerprint(x,2) for x in librarySet1]
nfps1 = len(fps1)
print "Number of molecules in library set: ", nfps1


#Similiarity searching

librarySet2 = [x for x in Chem.SDMolSupplier('/home/dsribar/Documents/DIO3/ROCS_2/hits_all.sdf')]
while librarySet2.count(None): librarySet2.remove(None)
fps2 = [GetMorganFingerprint(x,2) for x in librarySet2]
nfps2 = len(fps2)
print "Number of molecules in library set: ", nfps2


similiarity=[]
pairs = []
TanimotoCombo = []
    for x in range(0,len(fps2)):
        similiarity.append(DataStructs.TanimotoSimilarity(fps1[21],fps2[x]))
        pairs.append(str(librarySet1[21].GetProp("ID"))+" " + str(librarySet2[x].GetProp("ID")))
        TanimotoCombo.append(str(librarySet2[x].GetProp("ROCS_TanimotoCombo")))
d = {"Pairs" : pairs, "Similiarity" : similiarity, "TC":TanimotoCombo}
finalTable = pandas.DataFrame(d)

finalTable.to_csv("/home/dsribar/Documents/test.csv",index=False)

librarySet1[i].GetPropNames()[6]
