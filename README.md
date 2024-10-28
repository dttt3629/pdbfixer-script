# pdbfixer-script
A script to clean water and fix the non-stard residue with pdbfixer
## Prepare the protein for AI training from the original pdb file
##The script is used for clean water, replace non-stard res, add lacked res.
It comes to many issues so I put the issue and solve method below.
#Read them carefully to make sure they fit your job
# issue fix
## template = self.templates[residueName]  KeyError: 'UNK'
  pdbfixer.py line 511 add if residueName != 'UNK'

##KeyError: 'FOD' .Found in 6cxg
      i doubt it to be a variant res
      pdbfixer.py line 511 add if residueName not in ['UNK','FOD'] maybe is not quite right

##KeyError: 'NH2' .Found in 6wfz
      i doubt it to be a side of peptide linkage
      pdbfixer.py line 511 add if residueName not in ['UNK','FOD','NH2'] maybe is not quite right

## ValueError: could not convert string to float: ''
  del CRYSTAL1 line in pdb(cause it is empty,the author didnt fill the inf)   
 
## Cannot create a Context for a System with no particles
  a tough problem, found wrong in fixer.removeHeterogens(keepWater=False)
  cause the res is recognize wrong due to the strange form(exm TRYHH and the corect one is TRY H)
  use replace_from_i to fix pdbfile

## pdb like 7yar only keep the pos of CA, 
  use clean_pdb_ca to found and exclusive it
