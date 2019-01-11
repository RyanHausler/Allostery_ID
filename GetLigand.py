import scoria
import Bio
import numpy as np
from Bio.PDB import *
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from pypdb import *
import re
import scoria
from IPython.display import HTML
import os, sys
import copy

#Retrieves all ligands bound to a protein
#Deletes any ligands that are not within 5 A of query protein
"""
Given the filename of a molecule, creates a file containing ligands bound to the molecule

:params str filename: filename containing molecule

:returns: null
:rtype: null
"""
def get_ligand(filename, outputfile, query_file):
	#Load in query protein structure
	query = scoria.Molecule()
	query.load_pdb_into(query_file)
	#Load in moved protien
	full_protein = scoria.Molecule()
	full_protein.load_pdb_into(filename)
	#Gets a selection of all atoms that are protein or that won't be in a good ligand
	protein_selection = full_protein.select_atoms({'resname':["MSE", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", 
                "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT", "TIP", "DOD", "NAG", "EDO", "UNX", "SO4", 
                "ZN", "CA", "GOL", "CL", "NA", "MG", "ACT", "HEM", 
                "DC", "ZN", "FME", "DG", "DA", "CS", "DT", "CU", 
                "SO4", "PO4", "NA", "NH2", "CVM", "BMA", "FUL", "MES", 
                "GSH", "EPE", "F09", "IOD", 
                "PEG", "MAN", "FAD", "ACE", "MN", "NAP", "K", "UNK", 
                "EOH", "CD", "FUC", "FE"]})
	#invert the Selection so we have only ligands
	ligand_selection = full_protein.invert_selection(protein_selection)
	#Load selection into molecule
	ligand = full_protein.get_molecule_from_selection(ligand_selection)
	#Save pdb file
	ligand.save_pdb(outputfile)

	#Get each individual molecule
	#Store in list of lists (bound)
	bound = separate_into_molecules(filename)

	#For each molecule in the matched protein...
	mol_num = 0
	for mol in reversed(bound):
		#As long as it's not the actual protein
		if(mol_num != len(bound)):
		    #Convert mol to numpy array
		    mol_array = np.array(mol)
		    #create molecule from this selection
		    this_molecule = full_protein.get_molecule_from_selection(mol_array, serial_reindex=False)
			#If the molecule is further than 5 A from the protein...
		    if this_molecule.get_distance_to_another_molecules(query) > 5:     #NOT USING THE RIGHT PROTEIN
			    #remove mol from selection
			    for atom in reversed(mol):
				    full_protein.delete_atom(atom)
	mol = bound[0]
	for atom in reversed(mol):
		full_protein.delete_atom(atom)
	full_protein.save_pdb(outputfile)

#Get a list of individual ligands from a ligand selection
def separate_into_molecules(filename):
	#List of lists containing selections of each molecule
	bound = []
	#Read in each line of pdb file
	with open(filename) as myFile:
		content = myFile.readlines()
	mol = []
	atom_num = 0
	for line in content:
		#Check if next atom is hetatom
		ter = line[0:3]
		ter = str(ter)
		if(ter == 'END'):
			break
		if(ter == 'TER'):
			#Add the current molecule to bound and clear mol
			bound.append(mol)
			mol = []
			atom_num = atom_num - 1
		else:
			mol.append(atom_num)
		atom_num = atom_num + 1
	return(bound)

#Gets ligands from full protein
#These ligands will later be added to each chain separately
"""
Given the filename of a molecule, returns a selection of all
ligands in pdb file

:params str filename: pdb file containing protein and bound ligands

:returns numpy array of ligand selection 

"""
def get_ligand_as_strucutre(filename):
	protein = scoria.Molecule()
	protein.load_pdb_into(filename)
	atom_inf = protein.get_atom_information()
	#Get all atoms that are protein and selections we don't want
	protein_selection = protein.select_atoms({'resname':["MSE", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", 
                "SER", "THR", "TRP", "TYR", "VAL", "HOH", "WAT", "TIP", "DOD", "NAG", "EDO", "UNX", "SO4", 
                "ZN", "CA", "GOL", "CL", "NA", "MG", "ACT", "HEM", 
                "DC", "ZN", "FME", "DG", "DA", "CS", "DT", "CU", 
                "SO4", "PO4", "NA", "NH2", "CVM", "BMA", "FUL", "MES", 
                "GSH", "EPE", "F09", "IOD", 
                "PEG", "MAN", "FAD", "ACE", "MN", "NAP", "K", "UNK", 
                "EOH", "CD", "FUC", "FE"]})
	#invert the Selection
	ligand_selection = protein.invert_selection(protein_selection)
	#ligand = protein.get_molecule_from_selection(ligand_selection)
	return(ligand_selection)

#Gets selection of each chain with ligands bound
#Do before structural alignment
#Questions, talk to Patrick Ropp
"""
Given the filename of a molecule, returns a dictionary containing all 
chain molecules

:params str filename: filename containing molecule

:returns: dictionary of chain molecule selections
:rtype: dictionary
"""
def create_chain_molecules(pdbfile, ligand_selection):
    mol = scoria.Molecule()
    mol.load_pdb_into(pdbfile)
    #get_ligand_selections
    select_dict = mol.selections_of_chains()
    chain_mol = []
    for chain in select_dict.keys():
        selection = select_dict[chain]
		#append ligand selection to this
        selection = np.append(selection, ligand_selection)
        new_mol = mol.get_molecule_from_selection(selection)
        chain_mol.append(copy.deepcopy(new_mol))
    return(chain_mol)


#Gets the Amino Acid Sequence of each protein
#Any questions, ask Dr. Durrant
"""
Given the structure object of a molecule, obtains 
the sequence of amino acid

:params str filename: strucutre of protein

:returns: amino acid sequence
:rtype: string
"""
def get_seq(struc):
    ppb = PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(struc):
        sequence = sequence+pp.get_sequence()
        return sequence

#Function for aligning query protein with BLAST results
#Any questions, ask Dr. Durrant
"""
Given the structures the query protein and matched protein as well
as the name of the query protein, structurally aligns the matched protein
to the query protein and saves the to teh directory moved

:params structure template: strucutre of query protein 
        structure moving: strucutre of matched protein
        str query: name of query protein

:returns: null
:rtype: null
"""
def align_proteins(template, moving, query):
	cwd = os.getcwd()
	fastaFileName = cwd+'/data/Alignment_'+query+'/pdb'+moving+'.fasta'
	PDB_File_Name = cwd+'/data/Moved_'+query+'/pdb'+'moved'+moving+'.pdb'
	alnFileName = cwd+'/data/Alignment_'+query+'/pdb'+moving+'.aln'
	# Load in the moving PDB
	moving = PDBParser().get_structure(moving, cwd+'/data/Chains_'+query+'/pdb'+moving.lower()+".pdb")
	
	seq1 = get_seq(template)
	seq2 = get_seq(moving)
	
	# Make a fake fasta file
	header_line = ">protein1"
	header_line2 = ">protein2"
	t = open(fastaFileName,'w')
	t.write(header_line + "\n")
	t.write(str(seq1) + "\n")
	t.write("\n")
	t.write(header_line2 + "\n")
	t.write(str(seq2) + "\n")
	t.write("\n")
	t.close()

	# Use clustal to align the two sequences
	#info_path = os.path.dirname(os.path.abspath(__file__))
	#clustalw2_exe =  info_path + "/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
	clustalw2_exe = r"/usr/bin/clustalw"
	clustalw_cline = ClustalwCommandline(clustalw2_exe, infile=fastaFileName)
	stdout, stderr = clustalw_cline()
	
	# Load the alignment file into an alignment object
	alignment = AlignIO.read(alnFileName, "clustal")
	
	try:
		# Now map that alignment onto the two proteins
		# Tells us which residue in one structure is equivalent to the residue in the other structure
		# Based on the sequence alignment you performed with clustalw2.
		a = StructureAlignment(alignment, template, moving)
		maps = a.get_maps()
		# Now make lists, where entry x on each list corresponds to equivalent CA atoms
		ref_atoms = []
		moving_atoms = []
		for template_res in maps[0].keys():
			moving_res = maps[0][template_res]
			if not template_res is None and not moving_res is None:
				ref_atoms.append(template_res["CA"])
				moving_atoms.append(moving_res["CA"])
		# Now superimpose the structures by aligning those CAs
		super_imposer = Superimposer()
		super_imposer.set_atoms(ref_atoms, moving_atoms)
		super_imposer.apply(moving.get_atoms())

		# Save the aligned file
		io = PDBIO()
		io.set_structure(moving) 
		io.save(PDB_File_Name)
	except:
		print('Error. Too many gaps in the structure')