#!/usr/bin/env python3



import urllib
import pandas as pd
from io import StringIO
import requests


def get_uniprot_sequences(uniprot_ids):
	if type(uniprot_ids) == str:
		uniprot_ids = uniprot_ids.replace(" ", "").split(",")
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
	return df_fasta.assign(Query=df_fasta['Query'].str.split(',')).explode('Query')

if __name__ == "__main__":
	print (get_uniprot_sequences("Q29463"))