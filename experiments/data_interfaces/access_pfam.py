#!/usr/bin/env python3







import os
import numpy
import pickle
import xmltodict
import requests
from tqdm.autonotebook import tqdm

class pfam_records:
    def __init__(self, working_directory = "DEFAULT", prot_list = "EMPTY", use_cache_file=None):
        if working_directory== "DEFAULT":
            here, _ = os.path.split(os.path.abspath(os.path.realpath(__file__)))
            self.work_dir = here + "/cached_data"
            if not os.path.exists(self.work_dir):
                os.makedirs(self.work_dir)
        else:
            self.work_dir = os.path.abspath(working_directory)

        self.__pfam_record=dict()
        self.__proteins_with=dict()
        if prot_list != "EMPTY":
            self.add_prots(prot_list,use_cache_file)


    def get_prots_inset(self):
        return self.__pfam_record.keys()

    def get_domains_inset(self):
        return self.__proteins_with.keys()

    def get_proteins_with(self, dom):
        return self.__proteins_with[dom]

    def reduce_to_nonredunant(self, method = "segment30_string", prot_preference_score= None):
        if prot_preference_score ==None:
            s_set  = set()
            reduced_list = []
            for p in self.__pfam_record:
                s = pfam_records.__get_segment_string(self.__pfam_record[p])
                if s in s_set:
                    pass
                else:
                    s_set.add(s)
                    reduced_list.append(p)
            new_records_dict = dict()
            for p in reduced_list:
                new_records_dict[p] = self.__pfam_record[p]
            self.__pfam_record = new_records_dict

        else:
            representative_prot = dict()
            highest_score = dict()
            for p in self.__pfam_record:
                s = pfam_records.__get_segment_string(self.__pfam_record[p])
                if s not in highest_score:
                    highest_score[s] = prot_preference_score[p]
                    representative_prot[s] = p
                elif prot_preference_score[p] > highest_score[s]:
                    highest_score[s] = prot_preference_score[p]
                    representative_prot[s] = p
            reduced_list = list(representative_prot.values())
            new_records_dict = dict()
            for p in reduced_list:
                new_records_dict[p] = self.__pfam_record[p]
            self.__pfam_record = new_records_dict
        return reduced_list

    def get_pfam_record_for_prot(self, pid):
        try:
            ret =  self.__pfam_record[pid]
        except:
            raise Exception("not in this protein collection!")
        return ret



    def add_prots(self, pid_list,use_cache_file=None):
        if use_cache_file:
            # read the cache file which is a pickle of a dict like __pfam_record
            if os.path.exists(use_cache_file):
                with open(use_cache_file, "rb") as inf:
                    cached_record = pickle.load(inf)
            else:
                cached_record = dict()
        pftmp= self.work_dir + "/pfam_tmp"
        if type(pid_list) ==str:
            pid_list = [pid_list]
        if not os.path.exists(pftmp):
            os.makedirs(pftmp)
        for pid in tqdm(pid_list, desc="collecting domain annotations"):
            if pid in cached_record:
                p_dict = cached_record[pid]
            else:
                p_path = pftmp + "/" + pid + ".xml"
                if not os.path.exists(p_path):
                    url = "http://pfam.xfam.org/protein/" + pid + "?output=xml"
                    req = requests.get(url)
                    with open(p_path, "w") as outf:
                        outf.writelines(req.text)
                with open(p_path) as pf:
                    p_dict = xmltodict.parse(pf.read())
            self.__pfam_record[pid] = p_dict
            dlist = self.__extract_domain_list(p_dict)
            for d in dlist:
                if d in self.__proteins_with:
                    self.__proteins_with[d].add(pid)
                else:
                    self.__proteins_with[d] = {pid}

        if use_cache_file:
            cached_record.update(self.__pfam_record)
            with open(use_cache_file, "wb") as outf:
                pickle.dump(cached_record, outf)


    def get_domainlist_for_prot(self, pid):
        return pfam_records.__extract_domain_list(self.__pfam_record[pid])

    def get_segment_string(self,pid,blank_length_threshold=30):
        return pfam_records.__get_segment_string(self.__pfam_record[pid],blank_length_threshold)

    def get_singledom_prots(self, method = "segment_string30"):
        single_doms =[]
        for pid in self.__pfam_record.keys():
            pdict = self.__pfam_record[pid]
            if pfam_records.__check_single_domain(pdict, method=method):
                single_doms.append(pid)
        return single_doms

    def is_single_domain(self, pid,method="segment_string30"):
        return pfam_records.__check_single_domain(self.__pfam_record[pid],method=method)

    def extract_subset(self, pid_list):
        records = dict()
        pw= dict()
        for p in pid_list:
            r = self.__pfam_record[p]
            records[p] = r
            l = pfam_records.__extract_domain_list(r)
            for dom in l:
                if dom in pw:
                    pw[dom].add(p)
                else:
                    pw[dom] = {p}
        pfr = pfam_records()
        pfr.__pfam_record = records
        pfr.__proteins_with = pw
        return pfr

    #TODO: This should be modified it just does the job now but is not very well-defined
    def same_prots(self, prot_set): #this way of doing this function, assumes that we want to work with existing data in the two pfam_records object one of them is self an the other is the argument to the function.
        this_s_set = set()
        for p in self.__pfam_record:
            this_s_set.add(pfam_records.__get_segment_string(self.__pfam_record[p]))

        s_set = set()
        for p in prot_set.__pfam_record:
            s_set.add(pfam_records.__get_segment_string(prot_set.__pfam_record[p]))

        return this_s_set.intersection(s_set)


    @staticmethod
    def __get_segment_string(pdict, blank_length_threshold=30):
        if "pfam" not in pdict:
            return ""
        if "matches" not in pdict["pfam"]["entry"]:
            return ""
        p_seq = pdict["pfam"]["entry"]["sequence"]["#text"]
        p_len = len(p_seq)
        on_domain = [False] * p_len
        p_domains = pdict["pfam"]["entry"]["matches"]["match"]
        if type(p_domains) != list:
            p_domains = [p_domains]
        last_end = 0
        segments = []

        starts = [int(dom["location"]["@start"]) for dom in p_domains]
        idx = numpy.argsort(starts)
        sorted_array = numpy.array(p_domains)[idx]
        p_domains_sorted = list(sorted_array)
        for dom in p_domains_sorted:
            acc = dom["@accession"]
            begin = int(dom["location"]["@start"])
            end = int(dom["location"]["@end"])
            if begin - 1 - last_end > blank_length_threshold:
                segments += ["#", acc]
            else:
                segments += [acc]
            last_end = end
        if p_len - end > blank_length_threshold:
            segments += ["#"]
        return "_".join(segments)

    @staticmethod
    def __check_single_domain(pdict, method = "segment_string30"):
        seg_str = pfam_records.__get_segment_string(pdict,30)
        if seg_str=="":
            return False
        elif "_" in seg_str:
            return False
        elif "#" in seg_str:
            return False
        else:
            return True



    @staticmethod
    def __extract_domain_list(pdict):
        domain_list = []
        if "pfam" not in pdict:
            return []
        if "matches" in pdict["pfam"]["entry"]:
            p_domains = pdict["pfam"]["entry"]["matches"]["match"]
            if type(p_domains) != list:
                p_domains = [p_domains]
            for dom in p_domains:
                acc = dom["@accession"]
                domain_list.append(acc)
        return domain_list

if __name__ == "__main__":
    df = pfam_records()
    df.add_prots("P0A884")
    print(df.get_singledom_prots())

