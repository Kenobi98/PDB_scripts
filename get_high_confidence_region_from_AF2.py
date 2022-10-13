######### Obtain high confidence regions from AF2 model. #########
########## ------>"May the Force serve u well..." <------#########
##################################################################

############# One above all #############
##-------------------------------------##
import sys

import Bio
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import MMCIFParser

cif = sys.argv[1]         # CIF file path

models = MMCIFParser().get_structure('cif', cif)

high_confidence_regions = []
for model in models:
	for chain in model:
		for residue in chain.get_residues():
			for atom in residue:
				if atom.get_bfactor() > 80.0:
					high_confidence_regions.append(residue.id[1])				
				class AtomsSelect(Select):
				    def accept_atom(self, atom):
				        if atom.get_bfactor() > 80.0:
				            return 1
				        else:
				            return 0


io = PDBIO()
io.set_structure(models)
path = cif.split( ".cif" )[0]
io.save( f"{path}_confident.pdb", AtomsSelect() )

def ranges( positions ):
    nums = sorted( set( positions ) )
    gaps = [[s, e] for s, e in zip( positions, positions[1:] ) if s+1 < e]
    edges = iter(positions[:1] + sum(gaps, []) + positions[-1:])
    return list(zip(edges, edges))

confident_regions = ranges( high_confidence_regions )

print("\nRegions with a stretch of 2-10 residues: ")
[print(i) for i in confident_regions if i[1] - i[0] >1 and i[1] - i[0] <= 10]

print("\nRegions with a stretch of 11-20 residues: ")
[print(i) for i in confident_regions if (i[1] - i[0]) > 10 and (i[1] - i[0]) <= 20]

print("\nRegions with a stretch of 21-40 residues: ")
[print(i) for i in confident_regions if (i[1] - i[0]) > 20 and (i[1] - i[0]) <= 40]

print("\nRegions with a stretch of >40 residues: ")
[print(i) for i in confident_regions if (i[1] - i[0]) > 40]


	
		
