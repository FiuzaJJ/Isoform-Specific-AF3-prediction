from Bio.PDB import MMCIFParser
import os
import pandas as pd

def compute_contacts(path_to_mmcif, dist = 8 ):
    """Computes contacts from both chains named A and B of model 
    dist = 8 is the minimum distance to be counted as contact"""
    

    parser = MMCIFParser(QUIET=True) 
    structure = parser.get_structure("pdb_structure", path_to_mmcif)


    model = next(structure.get_models())  # first model
    chain1 = model["A"]
    chain2 = model["B"]
    contacts=[]    
    for res1 in chain1:
        for res2 in chain2:
            if res1.get_resname() == "GLY" or res2.get_resname() == "GLY":
                continue
                    
            else:
                if (res1["CB"]-res2["CB"])<dist:
                    contacts.append(((f"{res1.get_resname()}{res1.id[1]}",f"{res2.get_resname()}{res2.id[1]}"),res1["CB"]-res2["CB"]))
                    
    return contacts




def compute_af3_contacts(unzipped_dir,dist=8):

    """Computes contacts between the 2 strands of protein complexes on a unzipped folder provided from alphafold server
        returns for each key the contacts and amino acids involved"""
    results_dictionary={}
    
    # Walk through all subfolders and files
    for root, dirs, files in os.walk(unzipped_dir):
        
        for file in files:
            
            if "model_0" in file and file.endswith(".cif"):
                
                file_path = os.path.join(root, file)

                               
                name = file.removeprefix("fold_complex_")
                name = name.removesuffix("_model_0.cif")#clean up the names

                temp=compute_contacts(file_path,dist=dist)
                results_dictionary[name]=temp
                
 
    return results_dictionary


def contact_df(unzipped_dir,dist=8):

    """Computes contacts between the 2 strands of protein complexes on a unzipped folder provided from alphafold server
        returns for each key the contacts and amino acids involved
        """
    
    results=[]
    
    # Walk through all subfolders and files
    for root, dirs, files in os.walk(unzipped_dir):
        
        for file in files:
            
            if "model_0" in file and file.endswith(".cif"):
                
                file_path = os.path.join(root, file)

                               
                name = file.removeprefix("fold_complex_")
                name = name.removesuffix("_model_0.cif")#clean up the names

                temp=compute_contacts(file_path,dist=dist)
                results.append({
                    "complex_name": name,
                    "contact_count": len(temp)
                })
                
    df = pd.DataFrame(results)
    
    return df  


# def contact_df_v2(unzipped_dir,dist=8, SeqFile):

#     """Computes contacts between the 2 strands of protein complexes on a unzipped folder provided from alphafold server
#         returns for each key the contacts and amino acids involved
#         v2: Returns also collumn with contacts per aminoacid isoform length
#         collumn with contacts on diferential regions per isoform per aminoacid length"""
    
#     results=[]
#     seq_df=pd.read_csv(SeqFile)
    
#     # Walk through all subfolders and files
#     for root, dirs, files in os.walk(unzipped_dir):
        
#         for file in files:
            
#             if "model_0" in file and file.endswith(".cif"):
                
#                 file_path = os.path.join(root, file)

                               
#                 name = file.removeprefix("fold_complex_")
#                 name = name.removesuffix("_model_0.cif")#clean up the names

#                 temp=compute_contacts(file_path,dist=dist)
#                 results.append({
#                     "complex_name": name,
#                     "contact_count": len(temp)
#                 })
                
#     df = pd.DataFrame(results)
    
#     return df  


def translate(seq):
    """translates nucleotide sequence to protein"""
    print("running updated function")
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }

    protein = ""
    # clean up: remove whitespace/newlines, force uppercase
    seq = seq.upper().replace(" ", "").replace("\n", "")
    
    if len(seq) % 3 == 0:
        
        for i in range(0, len(seq), 3):

            codon = seq[i:i + 3].upper()
            
            protein += table.get(codon, "X")  # use X for unknown codons
    else:

        #### first block ignores last 1/2 elements of the codon to force a multiple by 3 solution
        # rest = len(seq) % 3 #here, substract the 1 or 2 last nucleotides to make the orf 3 dividable
        # seq=seq[:-rest]
        # for i in range(0, len(seq), 3):
        #     codon = seq[i:i + 3].upper()
                
        #     protein += table.get(codon, "X")  # use X for unknown codons


        ####
        protein+="NOT_3_DIVISIBLE"

    return protein