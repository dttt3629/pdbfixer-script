import os
from joblib import delayed,Parallel
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import pandas as pd
from tqdm import tqdm
import pdb

# issue fix
# template = self.templates[residueName]  KeyError: 'UNK'
#     pdbfixer.py line 511 add if residueName != 'UNK'

#KeyError: 'FOD' .Found in 6cxg, i doubt it to be a variant res
#      pdbfixer.py line 511 add if residueName not in ['UNK','FOD'] maybe is not quite right

#KeyError: 'NH2' .Found in 6wfz, i doubt it to be a side of peptide linkage
#      pdbfixer.py line 511 add if residueName not in ['UNK','FOD','NH2'] maybe is not quite right

# ValueError: could not convert string to float: ''
#     del CRYSTAL1 line in pdb(cause it is empty,the author didnt fill the inf)   
 
# Cannot create a Context for a System with no particles
#     a tough problem, found wrong in fixer.removeHeterogens(keepWater=False)
# cause the res is recognize wrong due to the strange form(exm TRYHH and the corect one is TRY H)
# use replace_from_i to fix pdbfile

# warning
# there are still many unslove problems in the script 
# example,7yar shows very strange, found it only made of CA,use clean_pdb_ca to found and exclusive it
# example,
def single_clean(pdbpath,save_str,keep_chain=None):
# pdbpath examp ./work/raw/1aru.pdb
# keep_chain list [A,B,C,D]
   try:
        fixer = PDBFixer(filename=pdbpath)
        #keep_main_chain
        if keep_chain is not None:
            delels = []
            for i1, i in enumerate(list(fixer.topology.chains())):
                if i.id not in keep_chain:
                    delels.append(i1)
            # for chainid in delels:
            fixer.removeChains(delels)

        # fix struct
        fixer.findMissingResidues()
        # only fix the missing res less than 20aa
        longres=[] 

        for i in fixer.missingResidues.keys():
             if len(fixer.missingResidues[i])>20:
                  longres.append(i)
        deled={}
        for i in longres:
             deled[i] = fixer.missingResidues[i].copy()
             del fixer.missingResidues[i]
        # dele long end res
        chainlength = []
        for i1, i in enumerate(list(fixer.topology.chains())):
              chainlength.append(len(list(i.residues())))
        endres = []
        for i in fixer.missingResidues.keys():
            if i[1] == 0 and len(fixer.missingResidues[i])>4:
                  endres.append(i)
            if i[1] == chainlength[i[0]]and len(fixer.missingResidues[i])>4:
                  endres.append(i)
        for i in endres:
             deled[i] = fixer.missingResidues[i].copy()
             del fixer.missingResidues[i]
        if endres ==[]:
            return 0    
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        #remove water
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        #pdb.set_trace()
        fixer.addMissingAtoms()

        fixer.addMissingHydrogens(7.0)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(save_str, 'w'),keepIds=True)
        # PDBFile.writeFile(fixer.topology, fixer.positions, open('out.pdb', 'w'),keepIds=True)
   except:
       0
def parallel_get(data,i):
        #try:
        
            Hchain = data.Hchain[i]
            Lchain = data.Lchain[i]
            pdb = data.pdb[i]
            atg = data.antigen_chain[i]
            if '|' in atg:
                atg = atg.replace(' ','')
                atg = atg.split('|')
            keep_chain = [Hchain,Lchain]
            keep_chain.extend(atg)
            path = './raw/'+pdb+'.pdb'
            savepath = './sab_clean/'+pdb
            for d in keep_chain:
                savepath+='_'
                savepath+=d
            # if not os.path.exists(savepath+'.pdb'):
            print(savepath)
            single_clean(path,save_str = savepath+'.pdb',keep_chain=keep_chain)
        # except:
        #     print(i)
def replace_res(path,chainid):
      f=open(path,'r')
      p = f.readlines()
      aalist = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 
                'TRP', 'SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN', 'ASP', 
                'GLU', 'LYS', 'ARG', 'HIS', 'MSE', 'CSO', 'PTR', 'TPO',
                'KCX', 'CSD', 'SEP', 'MLY', 'PCA', 'LLP']
      chainls = [i+i for i in chainid]
      p2 = []
      for i in p:
        for aa in aalist:
            for z in chainls:
                i= i.replace(aa+z,aa+' '+z[0])
        p2.append(i)
      f=open(path,'w')
      f.writelines(p2)
      f.close()
def replace_from_i(data,i):
            Hchain = data.Hchain[i]
            Lchain = data.Lchain[i]
            pdb = data.pdb[i]
            atg = data.antigen_chain[i]
            if '|' in atg:
                atg = atg.replace(' ','')
                atg = atg.split('|')
            keep_chain = [Hchain,Lchain]
            keep_chain.extend(atg)
            path = './raw/'+pdb+'.pdb'
            savepath = './sab_clean/'+pdb
            for d in keep_chain:
                savepath+='_'
                savepath+=d
            replace_res(path,keep_chain)
def clean_pdb_ca(data):
    wrongpdb = []
    for i in range(len(data)):
        pdb = data.pdb[i]
        path = './raw/'+pdb+'.pdb'
        fixer = PDBFixer(filename=path)
        a = list(fixer.topology.residues().atoms())
        if len(list(a[0].atoms()))<3:
               wrongpdb.append(path)
    return wrongpdb

if __name__ == '__main__':
    data = pd.read_csv('sab_sele.csv')
    # Parallel prepare, the openmm job is mp, so n_jobs is not sp to be huge 
    #Parallel(n_jobs = 2)(delayed(parallel_get)(i) for i in tqdm(range(len(data))))
    for i in tqdm(range(len(data))):
                   parallel_get(data,i)
    # a debug way to deal the wrong strange format pdb
    # while True:
    #     try:
    #         for i in tqdm(range(len(data))):
    #                parallel_get(data,i)
    #     # except OpenMMException:
    #     #     replace_from_i(i)
    #     #     # parallel_get(i)
    #     except Exception as e:
    #           print(e)
    #           replace_from_i(data,i)
    #           parallel_get(data,i)
    #           continue
    #         #   break
    #     break
        

