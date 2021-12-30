#!/usr/bin/env python3



import urllib
import pandas as pd
from io import StringIO
import requests
import os
import pickle

def get_uniprot_sequences(uniprot_ids, use_cache_file = "DEFAULT"):
	if use_cache_file == "DEFAULT":
		dir = os.path.dirname(os.path.abspath(__file__))
		use_cache_file = dir +"/cached_data/sequences_cache.pkl"

	if type(uniprot_ids) == str:
		uniprot_ids = uniprot_ids.replace(" ", "").split(",")
	cached_record = dict()
	if use_cache_file:
		# read the cache file which is a pickle of a dict like __pfam_record
		if os.path.exists(use_cache_file):
			with open(use_cache_file, "rb") as inf:
				cached_record = pickle.load(inf)
		else:
			dir,_  = os.path.split(use_cache_file)
			if not os.path.exists(dir):
				os.makedirs(dir)
	ids_to_fetch = list(set(uniprot_ids) - set(cached_record.keys()))
	fetched_seq = dict()
	if len(ids_to_fetch)>0:
		url = 'https://www.uniprot.org/uploadlists/'  # This is the webserver to retrieve the Uniprot data
		params = {
			'from': "ACC",
			'to': 'ACC',
			'format': 'tab',
			'query': " ".join(uniprot_ids),
			'columns': 'id,sequence'}
		data = urllib.parse.urlencode(params)
		data = data.encode('ascii')
		request = urllib.request.Request(url, data)
		with urllib.request.urlopen(request) as response:
			res = response.read()
		df_fasta = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
		df_fasta.columns = ["Entry", "Sequence", "Query"]
		# it might happen that 2 different ids for a single query id are returned, split these rows
		seq_df = df_fasta.assign(Query=df_fasta['Query'].str.split(',')).explode('Query')
		seq_df_pids = seq_df["Query"]
		seq_df_seqs = seq_df["Sequence"]
		for i in range(len(seq_df)):
			fetched_seq[seq_df_pids[i]] = seq_df_seqs[i]
	cached_record.update(fetched_seq)
	if use_cache_file:
		with open(use_cache_file, "wb") as outf:
			pickle.dump(cached_record, outf)
	return_dict = {x:cached_record[x] for x in uniprot_ids}
	return return_dict


if __name__ == "__main__":
	print (get_uniprot_sequences("Q29463"))