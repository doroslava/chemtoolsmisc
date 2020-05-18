from rdkit import Chem
import pandas
#Get the final list of selected compounds docked in CL097 and UMN
CL097SDFSet = [x for x in Chem.SDMolSupplier('/home/dsribar/Documents/tlr8/Repurposing/Docking/CL097/CL097_docked_repurposing_final.sdf')]
while CL097SDFSet.count(None): CL097SDFSet.remove(None)

UMNSDFSet = [x for x in Chem.SDMolSupplier('/home/dsribar/Documents/tlr8/Repurposing/Docking/UMN/UMN_docked_repurposing_final.sdf')]
while UMNSDFSet.count(None): UMNSDFSet.remove(None)



UMNdrugbankID = [x.GetProp("_Name") for x in UMNSDFSet]
CL097drugbankID = [x.GetProp("_Name") for x in CL097SDFSet]
UMNIds= []
CL097Ids=[]
for i in UMNdrugbankID:
    UMNIds.append(i.split("|")[0])

for i in CL097drugbankID:
    CL097Ids.append(i.split("|")[0])

AllIds = list(set(UMNIds)|set(CL097Ids))

drugbankTargets = pandas.read_csv("/home/dsribar/Documents/tlr8/Repurposing/drug_target.csv")
drugbankData = pandas.read_csv("/home/dsribar/Documents/tlr8/Repurposing/structure links.csv")


finalDataFrame = pandas.merge(drugbankTargets, drugbankData, on="DrugBank ID")

finalhits= finalDataFrame[finalDataFrame["DrugBank ID"].isin(AllIds)][["DrugBank ID","Name","Drug Groups","UniProt Name"]].drop_duplicates(subset="DrugBank ID")

finalhits.to_csv("/home/dsribar/Documents/tlr8/Repurposing/finalSelection.csv",index=False,sep="\t",header=True)