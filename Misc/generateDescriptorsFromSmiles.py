#Generates output file with calculated RDkit descriptors in .csv format from the input .csv file
#run with: python generateDescriptorsFromSmiles.py [input file] [output file]

import pandas
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from flask_mysqldb import MySQLdb


# SQL connection
db = MySQLdb.connect(host='mw3',
                             user='biasdb_user',
                             password='8asW3bWXLYNM26RM',
                             db='biasdb')
cur = db.cursor()
sql = "select * from ligand;"
cur.execute(sql)
result=cur.fetchall()


#inputLigandDataFrame = pandas.read_csv(sys.argv[1])

inputLigandDataFrame = pandas.read_csv("/home/dsribar/Documents/biasdb/biasdb_ligand_dataset_190807.csv")

outputLigandDataFrame = inputLigandDataFrame.copy()

#Generate descriptors and add them to corresponding columns

outputLigandDataFrame["Weight"] = [Chem.Descriptors.MolWt(Chem.MolFromSmiles(x)) for x in outputLigandDataFrame["smiles"]]
outputLigandDataFrame["LogP"]  = [Chem.Descriptors.MolLogP(Chem.MolFromSmiles(x)) for x in outputLigandDataFrame["smiles"]]
outputLigandDataFrame["NumHAcceptors"]=[Chem.Descriptors.NumHAcceptors(Chem.MolFromSmiles(x)) for x in outputLigandDataFrame["smiles"]]
outputLigandDataFrame["NumHDonors"]= [Chem.Descriptors.NumHDonors(Chem.MolFromSmiles(x)) for x in outputLigandDataFrame["smiles"]]
outputLigandDataFrame["NumRings"] =  [Chem.Descriptors.RingCount(Chem.MolFromSmiles(x)) for x in outputLigandDataFrame["smiles"]]
outputLigandDataFrame["TPSA"]= [Chem.Descriptors.TPSA(Chem.MolFromSmiles(x)) for x in outputLigandDataFrame["smiles"]]

#Check for duplicates or errors


#outputLigandDataFrame.to_sql("ligand",db,if_exists=="replace",index=False)
#outputLigandDataFrame.to_csv("/home/dsribar/Documents/biasdb/test_update.csv",index=False)

db.close()




