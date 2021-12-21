#!/usr/bin/env python3


import pandas as pd
from io import StringIO
import requests

max_per_request = 200
def get_SMILES(cid_list):
	num_batch = len(cid_list) // max_per_request
	df = pd.DataFrame()
	up_idx = 0    
	for i in range (num_batch):
		bot_idx = i * max_per_request
		up_idx = (i+1) * max_per_request
		# idxs  = list(range(bot_idx, up_idx))
		batch_list = cid_list[bot_idx:up_idx]
		url  = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/[[CIDLIST]]/property/CanonicalSMILES/CSV'.replace('[[CIDLIST]]',','.join(batch_list))
		result = requests.get(url)
		# temp_df = "ERROR"
		if result.ok:
			temp_df = pd.read_csv(StringIO(result.text), sep=",")
			df = pd.concat([df, temp_df], axis = 0) 
		else:
			print(url)
			print(result.text)
			raise Exception("error!")        
	# idxs  = list(range(up_idx, len(cid_list)))
	batch_list = cid_list[up_idx:]
	url  = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/[[CIDLIST]]/property/CanonicalSMILES/CSV'.replace('[[CIDLIST]]',','.join(batch_list))
	result = requests.get(url)
	if result.ok:
		temp_df = pd.read_csv(StringIO(result.text), sep=",")
		df = pd.concat([df, temp_df], axis = 0)
	else:
		print(url)
		print(result.text)

		raise Exception("error!")

	return df

if __name__=="__main__":
	print(get_SMILES(["448548"]))