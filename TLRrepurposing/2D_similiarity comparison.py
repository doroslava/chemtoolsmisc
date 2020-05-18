from os import listdir
from os.path import isfile
from os.path import join

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from e3fp.fingerprint.metrics import tanimoto, dice, cosine
from e3fp.fingerprint.fprint import Fingerprint, CountFingerprint
import numpy as np

import pandas as pd
#for histogram generation and manipulation
import seaborn
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#get the file names in the folder with fingerprints

onlyFiles = [f for f in listdir("/home/dsribar/Documents/tlr8/Repurposing/Fingerprints") if
             isfile(join("/home/dsribar/Documents/tlr8/Repurposing/Fingerprints", f))]
onlyNames = [s.strip('.fps.bz2') for s in onlyFiles]
#be sure to check if the name of the molecule actually corresponds to the right molecule. I didn´t finda a way to extract name information from
# Chem.SmilesMolSupplier, therefore the easiest was to rearrange the lines in drugs.smi file
suppl = Chem.SmilesMolSupplier('/home/dsribar/Documents/tlr8/Repurposing/drugs.smi',delimiter='\t',titleLine=False)

comparisonMatrix2D = []

for i in range(0, len(suppl)):
    row = []
    db1 = Fingerprint.from_rdkit(AllChem.GetMorganFingerprintAsBitVect(suppl[i],2,nBits=1024))
    for j in range(0, len(suppl)):
        db2 = Fingerprint.from_rdkit(AllChem.GetMorganFingerprintAsBitVect(suppl[j], 2,nBits=1024))
        row.append(round(tanimoto(db1, db2), 2))
    comparisonMatrix2D.append(row)
print(db2)
comparisonMatrix2D = pd.DataFrame(comparisonMatrix2D)
comparisonMatrix2D.columns = onlyNames
comparisonMatrix2D.index = onlyNames

seaborn.heatmap(comparisonMatrix2D, cmap="Blues", annot=True)
plt.yticks(rotation=0)
plt.xticks(rotation=60)
plt.title("Morgan fingerprints")
plt.show()







