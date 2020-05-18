from rdkit import Chem

librarySDFSet = [x for x in Chem.SDMolSupplier('/home/dsribar/Downloads/duplicates_check.sdf')]

libraryLength = len(librarySDFSet)
br =0

for i in range(0,libraryLength):
    for j in range(i+1,libraryLength):
        if (librarySDFSet[i].GetProp("Smiles")==librarySDFSet[j].GetProp("Smiles")):
            br =br+1;
            print(librarySDFSet[i].GetProp("Comp_ID"),librarySDFSet[j].GetProp("Comp_ID"))



from rdkit import Chem
a = []
suppl =  [x for x in Chem.SDMolSupplier('/home/dsribar/Downloads/duplicates_check.sdf')]
s = len(suppl)

for i in range(0, s):
    for j in range(i + 1, s):
        if (suppl[i].GetProp("canonical_smiles")==suppl[j].GetProp("canonical_smiles")):
            print(i,j)

suppl[2].GetPropNames()[3]