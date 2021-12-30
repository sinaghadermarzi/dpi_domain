#!/usr/bin/env python3
import pandas
import pandas as pd
from io import StringIO
import requests
import os
import pickle

max_per_request = 200

def get_SMILES_single_request(cid_list):
	url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/[[CIDLIST]]/property/CanonicalSMILES/CSV'.replace(
		'[[CIDLIST]]', ','.join(cid_list))
	result = requests.get(url)
	SMILES = dict()
	if result.ok:
		df = pd.read_csv(StringIO(result.text), sep=",")
		smiles_df_cids = df["CID"].values
		smiles_df_smiles = df["CanonicalSMILES"].values
		for i in range(len(smiles_df_cids)):
			SMILES[str(smiles_df_cids[i])] = smiles_df_smiles[i]
	else:
		print(url)
		print(result.text)
		raise Exception("error!")
	return SMILES


def get_SMILES(cid_list, use_cache_file = "DEFAULT"):
	ids_to_fetch= cid_list
	if use_cache_file == "DEFAULT":
		dir = os.path.dirname(os.path.abspath(__file__))
		use_cache_file = dir +"/cached_data/SMILES_cache.pkl"

	# read the cache file which is a pickle of a dict like __pfam_record
	if os.path.exists(use_cache_file):
		with open(use_cache_file, "rb") as inf:
			cached_record = pickle.load(inf)
	else:
		cached_record = dict()
		dir,_  = os.path.split(use_cache_file)
		if not os.path.exists(dir):
			os.makedirs(dir)
	ids_to_fetch= list(set(cid_list) - set(cached_record.keys()))

	if len(ids_to_fetch)>0:
		num_batch = len(ids_to_fetch) // max_per_request
		up_idx = 0
		for i in range (num_batch):
			bot_idx = i * max_per_request
			up_idx = (i+1) * max_per_request
			# idxs  = list(range(bot_idx, up_idx))
			batch_list = ids_to_fetch[bot_idx:up_idx]
			cached_record.update(get_SMILES_single_request(batch_list))
			with open(use_cache_file, "wb") as outf:
				pickle.dump(cached_record, outf)
		batch_list = ids_to_fetch[up_idx:]
		if len(batch_list)>0:
			cached_record.update(get_SMILES_single_request(batch_list))
			with open(use_cache_file, "wb") as outf:
				pickle.dump(cached_record, outf)

	return_dict = {x:cached_record[x] for x in cid_list}
	return return_dict

if __name__=="__main__":
	print(get_SMILES(["448548"]))