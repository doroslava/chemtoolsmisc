from rdkit.Chem import AllChem
from rdkit import Chem
import csv

def generate_combinatorial_library(reactant, reaction, dataset):
    
    """
    Function that generates virtual combinatorial libraries from the starting reactant and 
    library of building blocks and writes it to the current folder. Output file is in the
    SMILES format, where the second column is the ID of the respective building block.
    
    Arguments:
    reactant (char) : SMILES representation of the starting reactant
    reaction (char) : SMARTS pattern representation of the chemical reaction
    dataset  (char) : path to the file with building blocks, which are in the SMILES format
    
    Returns:
        None
    """
    rxn = AllChem.ReactionFromSmarts(reaction)
    
    with open(dataset) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            print ("Fragment is:", row[0])
            ps =rxn.RunReactants((Chem.MolFromSmiles(reactant),
                                  Chem.MolFromSmiles(row[0])))
            uniqps = {}
            for p in ps:
                smi = Chem.MolToSmiles(p[0])
                uniqps[smi] = p[0]
                print("Product is: ", sorted(uniqps.keys()))
            print("\n")
            for i in sorted(uniqps.keys()):
                with open("combinatorial_library_output.smi", "a+") as myfile:
                    myfile.write(i,)
                    myfile.write("\t")
                    myfile.write(row[1])
                    myfile.write("\n")
    