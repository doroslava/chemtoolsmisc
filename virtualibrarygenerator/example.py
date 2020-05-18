from librarygenerator import generate_combinatorial_library

#Specify reactants and reaction
reactant = 'CS(=O)(=O)C1=NC(=CC(=N1)C1=CC=CO1)C(F)(F)F'
reaction = '[SX4](=[OX1])(=[OX1])([#6:4][*:5])[#6].[#6:8][NX3;H2;!$(NC=[!#6]);!$(NC#[!#6]):9]>>[*:5]~[#6:4][N;H1:9][#6:8]'
building_blocks = "data/BB.smi"

#Generate virtual library and write it to the smiles file
generate_combinatorial_library(reactant,reaction,building_blocks)