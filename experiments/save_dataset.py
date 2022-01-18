#!/usr/bin/env python3
import os
import math
from data_interfaces.access_pubchem import get_SMILES
from data_interfaces.access_uniprot import get_uniprot_sequences
import hashlib



def hash_pairs(pair):
    pid =pair[0]
    cid =pair[1]
    pair_str =pid+cid
    return hashlib.md5(pair_str.encode(encoding='UTF-8')).hexdigest()

def save_dataset_csv(input_pairs,bindingdb_recs, output_path, only_with_ki=True):
    #this is to ensure a reproducible order of entries in the dataset while having shuffled them (prevent any hidden signal in the order of pairs)
    prot_cmp_pairs = sorted(input_pairs, key=hash_pairs)
    #prepare a csv file with inputs (SMILES strings for compounds and amino-acid sequences for proteins)
    #given a list of protein-compound pairs in which compounds are given by pubchem cid and proteins are given by uniprot id
    cid_list = set()
    pid_list = set()
    for (p , c) in prot_cmp_pairs:
        cid_list.add(c)
        pid_list.add(p)

    SMILES = get_SMILES(cid_list)
    seq = get_uniprot_sequences(pid_list)
    with open (output_path, "w") as outf:
        outf.writelines("protein_uniprot_id,compound_pubchem_cid,protein_sequence,compound_SMILES,binary_interaction,affinity[-log10(Ki/10e9)]\n")
        for (prot, compound) in prot_cmp_pairs:
            aff = bindingdb_recs.get_pAff_avg(prot, compound, "Ki")
            if type(aff)==str:
                aff_str = "none"
            else:
                aff_str= str(bindingdb_recs.get_pAff_avg(prot, compound,"Ki"))
                
            if (aff_str=="none") and only_with_ki:
                pass
            else:
                bin_int = str(int(bindingdb_recs.get_bin_interaction((prot,compound))))
                line = prot +"," +compound+ ","+seq[prot]+ ","+ SMILES[compound] + ","+bin_int+","+aff_str+"\n"
                outf.writelines(line)




from data_interfaces.access_bindingdb import bindingdb_records


if __name__== "__main__":
    bdb = bindingdb_records()
    bdb.read_from_pickle("./data_interfaces/cached_data/bindingdb2021m11_bin_Ki_1uM_30uM.pkl")
    pairs = [("P00742","6102595"),
             ("P09871",	"44333833"),
             ("P00740",	"11684611")]
    save_dataset_csv(pairs, bdb, "sina.csv")
    a= 7