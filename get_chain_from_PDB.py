######### Obtain a specific chain from PDB file. #########
###### ------>"May the Force serve u well..." <------#####
##########################################################

############# One above all #############
##-------------------------------------##
import sys

import Bio
from Bio.PDB import PDBParser, PDBIO, Select


##Input files ###########################
##-------------------------------------##
pdb = sys.argv[1]         # PDB file path
chain_ = sys.argv[2]         # Chain to be extracted.

models = PDBParser().get_structure('pdb', pdb)

high_confidence_regions = []
for model in models:
	for chain in model:
		class ChainSelect(Select):
		    def accept_chain(self, chain):
		        if chain.id == chain_:
		            return 1
		        else:
		            return 0

io = PDBIO()
io.set_structure(models)
path = pdb.split( ".pdb" )[0]
io.save( f"{path}_{chain_}.pdb", ChainSelect() )
