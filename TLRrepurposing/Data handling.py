#Script for modifying the data downloaded from Drugbank
import pandas

#Takes the input dataframe, where every target for specific molecule is in different row, and paste them in on one UniProt
# Name field- one row-one molecule (DrugBank ID)
drugbankDataFrame = pandas.read_csv("/home/dsribar/Downloads/uniprot links.csv")
drugbankDataFrameModified = drugbankDataFrame.groupby('DrugBank ID')['UniProt Name'].apply(', '.join).reset_index()
drugbankDataFrameModified.to_csv("/home/dsribar/Documents/tlr8/Repurposing/drug_target.csv",index=False)

