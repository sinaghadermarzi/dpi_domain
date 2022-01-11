from evaluate import *
import numpy
from tqdm.autonotebook import tqdm
def convert_aff_to_p(aff):
	if aff==0:
		outval = "invalid"
	else:
		outval= 9-numpy.log10(aff)
	return outval

def baseline_results(bdb,test_set_df,target):
	#TODO: clean up this (less repetition)
	df = test_set_df
	valid_aff = df["affinity[-log10(Ki/10e9)]"]!="none"
	valid_bin = df["binary_interaction"]!="none"
	aff_col = df.loc[valid_aff,"affinity[-log10(Ki/10e9)]"].astype(float).values
	bina = df.loc[valid_bin,"binary_interaction"].astype(int).values
	if target=="affinity":		
		prot_col = df.loc[valid_aff,"protein_uniprot_id"].values
		compound_col = df.loc[valid_aff,"compound_pubchem_cid"].values
	elif target=="binary":		
		prot_col = df.loc[valid_bin,"protein_uniprot_id"].values
		compound_col = df.loc[valid_bin,"compound_pubchem_cid"].values
	
	avg_aff = numpy.mean(aff_col)
	if target=="affinity":
		BP_CIs = []
		BP_SPs = []
		BP_MSEs = []
		BP_MAEs = []
		BP_PEAs = []
		BD_CIs = []
		BD_SPs = []
		BD_MSEs = []
		BD_MAEs = []
		BD_PEAs = []
	elif target=="binary":
		BP_AUCs = []
		BP_AUPRs = []
		BP_MCCs = []
		BP_F1 = []
		BP_SENs = []
		BP_SPECs = []
		BD_AUCs = []
		BD_AUPRs = []
		BD_MCCs = []
		BD_F1 = []
		BD_SENs = []
		BD_SPECs = []


	# prots = set()
	# compounds = set()
	# for i in range(len(prot)):
	# 	prots.add(prot[i])
	# 	compounds.add(drug[i])

	# cmp_aff_list = dict()
	# prot_aff_list = dict()
	# # for (p,c) in bdb.get_all_bin_pairs():
	# # 	if p in prots:
	# # 		aff_list= bdb.get_affinities(p,c)["Ki"]
	# # 		if p not in prot_aff_list:
	# # 			prot_aff_list[p] = []
	# # 		prot_aff_list[p] += aff_list
	# # 		cmp_aff_list +=  dict()

	# # 	if c in compounds:
	# # 		aff_list= bdb.get_affinities(p,c)["Ki"]
	# # 		if c not in cmp_aff_list:
	# # 			cmp_aff_list[c] = []
	# # 		cmp_aff_list[c] += aff_list
	aff_pool_p = dict()
	aff_pool_c = dict()
	for i in range(len (prot_col)):
		p = prot_col[i]
		c = compound_col[i]
		pairs = bdb.get_pairs_for_drug(c)
		other_pairs = [(pp,dd) for (pp,dd) in pairs if pp != p]
		
		all_pki =[]
		for (pp,dd) in other_pairs:
			all_pki += [convert_aff_to_p(x) for (x,y) in bdb.get_affinities(pp,dd)["Ki"]]
		all_pki = list(filter(("invalid").__ne__, all_pki))
		aff_pool_c[(p,c)] = all_pki

		pairs = bdb.get_pairs_for_prot(p)
		other_pairs = [(pp,dd) for (pp,dd) in pairs if dd != c]

		all_pki =[]
		for (pp,dd) in other_pairs:
			all_pki += [convert_aff_to_p(x) for (x,y) in bdb.get_affinities(pp,dd)["Ki"]]
		all_pki = list(filter(("invalid").__ne__, all_pki))
		aff_pool_p[(p,c)] = all_pki		


	for j in tqdm(range(100),desc="repetition"):
		baseline3_C_col = [None] * len (prot_col)
		baseline3_P_col = [None] * len (prot_col)
		for i in range(len (prot_col)):
			p = prot_col[i]
			c = compound_col[i]

			all_pki = aff_pool_c[(p,c)]
			if len(all_pki)==0:
				sel_pki = avg_aff
			else:
				sel_pki = numpy.random.choice(all_pki,1)[0]
			baseline3_C_col[i] = sel_pki
			#for baseline3_P, we take all the Kis of the given protein except those that are with the given drug
			all_pki = aff_pool_p[(p,c)]
			if len(all_pki)==0:
				sel_pki = avg_aff
			else:
				sel_pki = numpy.random.choice(all_pki,1)[0]
			baseline3_P_col[i] = sel_pki
		try:			
			if target=="binary":
				auc,aupr,mcc,f1, sensitivity, specificity= report_for_binary(bina,baseline3_P_col,print_report=False)
				BP_AUCs.append(auc)
				BP_AUPRs.append(aupr)
				BP_MCCs.append(mcc)
				BP_F1.append(f1)
				BP_SENs.append(sensitivity)
				BP_SPECs.append(specificity)
				auc,aupr,mcc,f1, sensitivity, specificity = report_for_binary(bina,baseline3_C_col,print_report=False)
				BD_AUCs.append(auc)
				BD_AUPRs.append(aupr)
				BD_MCCs.append(mcc)
				BD_F1.append(f1)
				BD_SENs.append(sensitivity)
				BD_SPECs.append(specificity)
			elif target=="affinity":
				ci,mae,mse,pearson,spearman = report_for_affinity(aff_col,baseline3_P_col,print_report=False)
				BP_CIs.append(ci)
				BP_SPs.append(spearman)
				BP_MSEs.append(mse)
				BP_MAEs.append(mae)
				BP_PEAs.append(pearson)
				ci,mae,mse,pearson,spearman  = report_for_affinity(aff_col,baseline3_C_col,print_report=False)
				BD_CIs.append(ci)
				BD_SPs.append(spearman)
				BD_MSEs.append(mse)
				BD_MAEs.append(mae)
				BD_PEAs.append(pearson)
		except:
			pass

	if target=="binary":
		print("Binary:\n-------------------------")
		print("results for the baseline3_P")
		print("AUC:\t",concise(BP_AUCs))
		print("AUPR:\t",concise(BP_AUPRs))
		print("MCC:\t",concise(BP_MCCs))
		print("F1:\t",concise(BP_F1))
		print("SEN:\t",concise(BP_SENs))
		print("SPEC:\t",concise(BP_SPECs))
		print("\nresults for the baseline3_C")
		print("AUC:\t",concise(BD_AUCs))
		print("AUPR:\t",concise(BD_AUPRs))
		print("MCC:\t",concise(BD_MCCs))
		print("F1:\t",concise(BD_F1))
		print("SEN:\t",concise(BD_SENs))
		print("SPEC:\t",concise(BD_SPECs))
	elif target=="affinity":
		print("\nAffinity:")
		print("results for the baseline3_P")
		print("CI:\t",concise(BP_CIs))
		print("MAE:\t",concise(BP_MAEs))
		print("MSE:\t",concise(BP_MSEs))
		print("SP:\t",concise(BP_SPs))
		print("PEA:\t",concise(BP_PEAs))
		print("\nresults for the baseline3_C")
		print("CI:\t",concise(BD_CIs))
		print("MAE:\t",concise(BD_MAEs))
		print("MSE:\t",concise(BD_MSEs))
		print("SP:\t",concise(BD_SPs))
		print("PEA:\t",concise(BD_PEAs))
		

