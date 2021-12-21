#!/usr/bin/env python3




import pickle
from tqdm.dask import TqdmCallback
from tqdm.autonotebook import tqdm
import math

import dask.dataframe as dd

import os



# this is a class to encapsulate protein-compound interaction data in a bindingdb database.

class bindingdb_records:
    def __init__(self):
        self.__records = None
        self.__prots_by_drug = None
        self.__drugs_by_protein = None

        self.__bin_interaction = None
        self.__bin_ligandsof = None
        self.__bin_targetsof = None


###### operations with the database (read/write/binarize/...)
    def __extract_records_from_dataframe(self, df):
        self.__records = dict()
        self.__prots_by_drug = dict()
        self.__drugs_by_protein = dict()
        # df_dict = df.to_dict("records")
        ## TODO: there should be a more efficient way of doing this loop. Looks like iterrows is very inefficient
        # for row in tqdm(df_dict, total=len(df), desc="extracting records"):
        for index, row in tqdm(df.iterrows(), total=len(df), desc="extracting records"):
            pid = row["UniProt (SwissProt) Primary ID of Target Chain"]
            valid_pid = (type(pid)==str) and (pid != "") and ("," not in pid)
            did = row["PubChem CID"]
            valid_did = (type(did)==str) and (did != "") and ("," not in did)
            single_chain = row["Number of Protein Chains in Target (>1 implies a multichain complex)"]=="1"
            valid_pair = (valid_did and valid_pid and single_chain)
            if valid_pair:
                record = dict()
                record["BindingDB Reactant_set_id"] = row["BindingDB Reactant_set_id"]
                record["Ki"] = row["Ki (nM)"]
                record["Kd"] = row["IC50 (nM)"]
                record["IC50"] = row["Kd (nM)"]
                record["EC50"] = row["EC50 (nM)"]
                record["kon (M-1-s-1)"] = row["kon (M-1-s-1)"]
                record["koff (s-1)"] = row["koff (s-1)"]

                if (pid, did) in self.__records:
                    self.__records[(pid, did)].append(record)
                else:
                    self.__records[(pid, did)] = [record]

                if pid in self.__drugs_by_protein:
                    self.__drugs_by_protein[pid].add(did)
                else:
                    self.__drugs_by_protein[pid]= {did}

                if did in self.__prots_by_drug:
                    self.__prots_by_drug[did].add(pid)
                else:
                    self.__prots_by_drug[did]= {pid}
        return None


    def read_from_csv(self, csv_path, sepa ='\t', df_pickle_dir = "./cached_data/"):
        # reads the dataframe from a csv and writes it down to a pickle file
        fname = ".".join(csv_path.split("/")[-1].split(".")[:-1])
        cb = TqdmCallback(desc="reading csv")
        cb.register()
        df = dd.read_csv(csv_path,blocksize=100e6, sep=sepa,dtype="str",on_bad_lines="skip").compute()
        cb.unregister()
        ## for now this does nothing, we don't need this 
        # df_pickle_dir = None
        # if df_pickle_dir != None:
        #     if df_pickle_dir == "./cached_data/":
        #         if not os.path.exists("./cached_data/"):
        #             os.makedirs("./cached_data/")
        #     df_pickle_path = df_pickle_dir+ fname+"_dataframe.pickle"
        #     with open (df_pickle_path, "wb") as outf:
        #         pickle.dump(df, outf)
        self.__extract_records_from_dataframe(df)


    ## for now this does nothing, we don't need this 
    # def read_from_dataframe_pickle(self,df_pickle_path):
    #     with open (df_pickle_path, "rb") as inf:
    #         df= pickle.load(inf)
    #     self.__extract_records_from_dataframe(df)
    #     return None

    def read_from_pickle(self,pickle_path):
        with open (pickle_path, "rb") as inf:
            [self.__records, self.__prots_by_drug, self.__drugs_by_protein, self.__bin_interaction,self.__bin_ligandsof,self.__bin_targetsof] = pickle.load(inf)
        return None

    def save_to_pickle(self,pickle_path):
        with open (pickle_path, "wb") as outf:
            pickle.dump([self.__records, self.__prots_by_drug, self.__drugs_by_protein, self.__bin_interaction,self.__bin_ligandsof,self.__bin_targetsof], outf)
        return None


    def binarize_db(self, l_threshold=1000, h_threshold=30000, use_measure = "Ki"):
        self.__bin_interaction = dict()
        self.__bin_ligandsof = dict()
        self.__bin_targetsof = dict()
        if use_measure=="use_all4":
            pass # to be implemented later if necessary for now we focus on Ki
        else:
            idx = 0
            for (p, d) in tqdm(self.__records, desc= "binarizing affinities"):
                idx += 1
                recs = self.__records[(p,d)]
                res = bindingdb_records.__binarize_use_specific(pair_records=recs,measure_name=use_measure,l_threshold=l_threshold,h_threshold=h_threshold)
                if type(res) == bool:
                    self.__bin_interaction[(p, d)] = res
                    if res:
                        if p in self.__bin_ligandsof:
                            self.__bin_ligandsof[p].add(d)
                        else:
                            self.__bin_ligandsof[p] = {d}

                        if d in self.__bin_targetsof:
                            self.__bin_targetsof[d].add(p)
                        else:
                            self.__bin_targetsof[d] = {p}
        return None



    ##### query operations with database records
    def get_pairs_for_drug(self, cid):
        return [(p, cid) for p in self.__prots_by_drug[cid]]

    def get_pairs_for_prot(self, pid):
        return [(pid, d) for d in self.__drugs_by_protein[pid]]

    def get_records_for_pair(self, pd_pair):
        return self.__records[pd_pair]

    def get_all_pairs(self):
        return list(self.__records.keys())

    def get_all_bin_pairs(self):
        return list(self.__bin_interaction.keys())

    def get_all_prots(self):
        return list(self.__drugs_by_protein.keys())

    def get_all_drugs(self):
        return list(self.__prots_by_drug.keys())

##### query operations with binary interaction
    def get_ligands_of(self,pid):
        if pid in self.__bin_ligandsof:
            return self.__bin_ligandsof[pid]
        else:
            return set()

    def get_targets_of(self,cid):
        if cid in self.__bin_ligandsof:
            return self.__bin_ligandsof[cid]
        else:
            return set()

    def get_bin_interaction(self,pd_pair):
        if pd_pair in self.__bin_interaction:
            return self.__bin_interaction[pd_pair]
        else:
            return "pair not found in db!"


##### query operations with affinities
    def get_affinities(self,p,d):
        outdict = dict()
        records = self.__records[(p, d)]
        for aff_type in ["Ki", "Kd", "IC50", "EC50"]:
            measurements = []
            for rec in records:
                lower_than = False
                higher_than = False
                field = rec[aff_type]
                if type(field)==str:
                    field = field.strip()
                    if ">" in field:
                        field = field.replace(">","")
                        higher_than = True
                    if "<" in field:
                        field = field.replace("<","")
                        higher_than = True
                if (higher_than and lower_than):
                    pass
                elif type(field)==float:
                    if math.isnan(field):
                        pass
                else:
                    try:
                        num=float(field)
                        if higher_than:
                            measurements.append((num,">"))
                        elif lower_than:
                            measurements.append((num,"<"))
                        else:
                            measurements.append((num,"="))
                    except:
                        pass
            outdict[aff_type] =measurements
        return outdict


    def get_pAff_avg(self,pid,cid,aff_type):
        aff_dict = self.get_affinities(pid, cid)
        vals = [ v for (v,rel) in aff_dict[aff_type]]
        aff= -math.log10((sum(vals)/len(vals))/10e9)
        return aff


### Utility functions
    @staticmethod
    def __binary_from_str(x, l_threshold, h_threshold):
        try:
            if type(x) != str:
                return "nonstring"
            else:
                x_c = x.replace(">", "")
                x_c = x_c.replace("<", "")
                x_c = x_c.replace(" ", "")
                val = float(x_c)
                if val == 0:
                    return "zero"
                if val > h_threshold:
                    val_bin = False
                elif val < l_threshold:
                    val_bin = True
                else:
                    return "middle"
                if val_bin and (">" in x):
                    return "invalid"
                if (not val_bin) and ("<" in x):
                    return "invalid"
                return val_bin
        except:
            return "error"


    @staticmethod
    def __binarize_use_several_any(pair_records, l_threshold, h_threshold):
        # to be implemented if neccessary
        # num_pos = 0
        # num_neg = 0
        # for i in range(len(pair_records)):
        #     bin_Ki = bindingdb_records.__binary_from_str(pair_records[i]["Ki"], l_threshold, h_threshold)
        #     bin_IC50 = bindingdb_records.__binary_from_str(pair_records[i]["IC50"], l_threshold, h_threshold)
        #     if (type(bin_Ki) == bool) and (type(bin_IC50) == bool):
        #         if bin_Ki and bin_IC50:
        #             num_pos += 1
        #         elif (not bin_Ki) and (not bin_IC50):
        #             num_neg += 1
        #     elif type(bin_Ki) == bool:
        #         if bin_Ki:
        #             num_pos += 1
        #         else:
        #             num_neg += 1
        #     elif type(bin_IC50) == bool:
        #         if bin_IC50:
        #             num_pos += 1
        #         else:
        #             num_neg += 1
        # num_all = num_pos + num_neg
        # if num_all == 0:
        #     return "undecided_none"
        # if num_neg == 0:
        #     return True
        # if num_pos == 0:
        #     return False
        return "undecided_conflict"


    @staticmethod
    def __binarize_use_several_all(pair_records, l_threshold, h_threshold):
        #to be implemented if necessary
        # num_pos = 0
        # num_neg = 0
        # for i in range(len(pair_records)):
        #     bin_Ki = bindingdb_records.__binary_from_str(pair_records[i]["Ki"], l_threshold, h_threshold)
        #     bin_IC50 = bindingdb_records.__binary_from_str(pair_records[i]["IC50"], l_threshold, h_threshold)
        #     if (type(bin_Ki) == bool) and (type(bin_IC50) == bool):
        #         if bin_Ki and bin_IC50:
        #             num_pos += 1
        #         elif (not bin_Ki) and (not bin_IC50):
        #             num_neg += 1
        #     elif type(bin_Ki) == bool:
        #         if bin_Ki:
        #             num_pos += 1
        #         else:
        #             num_neg += 1
        #     elif type(bin_IC50) == bool:
        #         if bin_IC50:
        #             num_pos += 1
        #         else:
        #             num_neg += 1
        # num_all = num_pos + num_neg
        # if num_all == 0:
        #     return "undecided_none"
        # if num_neg == 0:
        #     return True
        # if num_pos == 0:
        #     return False
        return "undecided_conflict"

    @staticmethod
    def __binarize_use_specific(pair_records, measure_name, l_threshold, h_threshold):
        #this function takes the id of a drug and a protein and based on the records of affinity measurements in the database, returns a binary interaction (True/False) or an error string
        num_pos = 0
        num_neg = 0
        for i in range(len(pair_records)):
            bin_val = bindingdb_records.__binary_from_str(pair_records[i][measure_name], l_threshold, h_threshold)
            if (type(bin_val) == bool):
                if bin_val:
                    num_pos += 1
                else:
                    num_neg += 1
        num_all = num_pos + num_neg
        if num_all == 0:
            return "undecided_none"
        elif num_neg == 0:
            return True
        elif num_pos == 0:
            return False
        else:
            return "undecided_conflict"








if __name__=="__main__":
    db = bindingdb_records()
    if (os.path.exists("./cached_data/BindingDB_All_2021m10.tsv_records.pickle")):
        db.read_from_pickle("./cached_data/BindingDB_All_2021m10.tsv_records.pickle")
    else:
        db.read_from_csv("../data/bindingdb/db2021m0.tsv")
        db.binarize_db(l_threshold=1000, h_threshold=30000)
        db.save_to_pickle("./cached_data/BindingDB_All_2021m10.tsv_records.pickle")