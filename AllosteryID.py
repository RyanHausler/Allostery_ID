import Bio
from Bio.PDB import *
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from pypdb import *
import re
import scoria
from IPython.display import HTML
import os, sys
import pprint
import glob
from GetLigand import get_ligand, align_proteins, get_seq, create_chain_molecules, get_ligand_as_strucutre
import copy
#-- AllosteryID
#-- Ryan Hausler
#Input: PDB ID of query protein
#Output: Aligned/unaligned PDB structures of similar BLAST results
#	 Files containing only the ligands that were bound to each BLAST result

#Create a pdb list
cwd = os.getcwd()
pdbl = PDBList()

#Ask for query protien
protein=input("Please enter query protein: ")
query=protein
proteinName=cwd+'/data/Original_'+query+'/pdb'+protein+'.ent'

#Make directory for Alignment files
if not os.path.exists(cwd+'/data/Alignment_'+protein+'/'):
    os.makedirs(cwd+'/data/Alignment_'+protein+'/')
#Make directory for moved pdb files
if not os.path.exists(cwd+'/data/Moved_'+query+'/'):
    os.makedirs(cwd+'/data/Moved_'+query+'/')
#Make directory for original pdb files
if not os.path.exists(cwd+'/data/Original_'+query+'/'):
    os.makedirs(cwd+'/data/Original_'+query+'/')
#Make Directory for Chains
if not os.path.exists(cwd+'/data/Chains_'+query+'/'):
    os.makedirs(cwd+'/data/Chains_'+query+'/')
#Make Directory for Ligands
if not os.path.exists(cwd+'/data/Ligands_'+query+'/'):
    os.makedirs(cwd+'/data/Ligands_'+query+'/')
#Create parser object
parser = PDBParser()

#Obtain query protein ID
fir=cwd+'/data/Original_'+query+'/'
pdbl.retrieve_pdb_file(protein, file_format='pdb', pdir=fir)

#Download structure of query protein
template = parser.get_structure(protein, proteinName)

#Runs BLASTP on pdb id
#Uses pypdb package
blast_results = get_blast2(protein, chain_id='A')

#Parses through BLAST results
#Obtains matched protein ID and score
#Saves each tuple into a list
#example [(id, score), (id, score),...]
myList=[]
for protein in range(len(blast_results[0])):
	results = str(blast_results[1][protein])
	id_start = results.find('</a>')+4
	id_end = results.find(':', 0)
	score_start = results.find('Score = ')+8
	score_end = results.find(' bits', score_start)
	tup = (results[id_start:id_end], int(float(str(results[score_start:score_end]))))
	myList.append(tup)

#Download all BLASTp matches
for match in myList:
    protein = match[0]
    pdbl.retrieve_pdb_file(protein, file_format='pdb', pdir=fir)

#Align only the BLAST results with a score above 250 to query protein
#Downloads the aligned/moved proteins
#Uses get_seq and align_proteins methods in GetLigands.py
for match in myList:
    if match[1] > 250:
        #Name of ith protein on list
        compare_protein = match[0]
        #Get the structure of the ith protein on the list
        moving = cwd+'/data/Original_'+query+'/pdb'+compare_protein.lower()+".ent"
        #Get all ligands from molecule
        ligand_selection = get_ligand_as_strucutre(moving)
        print(moving)
        print(ligand_selection)
        #Create a chain molecule from this structure
        chain_mol = create_chain_molecules(moving, ligand_selection)
        i = 1
        #Save each chain and align with the query protein
        for chain in chain_mol:
            name = compare_protein.lower()+str(i)
            save_name = cwd+'/data/Chains_'+query+'/pdb'+name+'.pdb'
            chain.save_pdb(save_name)
            #Align each chain with template
            align_proteins(template, name, query)
            i=i+1

#For every moved chain, get all of the ligands and 
#place them in individual PDB files
#Uses get_ligand methods in GetLigands.py
for filename in glob.glob(cwd+'/data/Moved_'+query+'/'+"*1.pdb"):
    filename2 = filename.lstrip(cwd+'/data/Moved_'+query+'/')
    outputfilename = cwd+'/data/Ligands_'+query+'/'
    get_ligand(filename, outputfilename+'ligand_'+filename2, proteinName)   