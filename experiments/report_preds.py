import dask.dataframe as dd
from evaluate import report_for_affinity,report_for_binary
import numpy





def evaluate_pred(dataset_path,pred_file, report_affinty= True):
    df = dd.read_csv(dataset_path,dtype="str",on_bad_lines="skip").compute()
    pred_dict = dict()

    with open(pred_file) as f:
        for line in f:
            line = line.strip()
            line = line.split(",")
            if len(line)<2:
                continue
            else:
                pred_dict[line[0]] = float(line[1])   
    avg_pred = numpy.mean(list(pred_dict.values()))


    ## for binary
    print("report for binary:\n")
    valid_row = (df["binary_interaction"] != "none")
    prot_col = df.loc[valid_row, "protein_uniprot_id"].values
    cmp_col = df.loc[valid_row, "compound_pubchem_cid"].values
    label_col = df.loc[valid_row, "binary_interaction"].astype(int)
    
    
    pred_col = [None]*len(prot_col)
    
    for i in range(len(prot_col)):
        pair_str = prot_col[i]+"_"+cmp_col[i]
        if pair_str in pred_dict:
            pred_col[i] = pred_dict[pair_str]
        else:
            pred_col[i] = avg_pred
            print(pair_str+" NOT FOUND!")
    print("\n",end="")
    report_for_binary(label_col,pred_col)
    print("-----------------------\n")
    #for affinity
    if report_affinty:
        print("report for affinity:\n")
        valid_row = (df["affinity[-log10(Ki/10e9)]"] != "none")
        prot_col = df.loc[valid_row, "protein_uniprot_id"].values
        cmp_col = df.loc[valid_row, "compound_pubchem_cid"].values
        aff_col = df.loc[valid_row,"affinity[-log10(Ki/10e9)]"].astype(float)
        pred_col = [None]*len(prot_col)
    
        for i in range(len(prot_col)):
            pair_str = prot_col[i]+"_"+cmp_col[i]
            if pair_str in pred_dict:
                pred_col[i] = pred_dict[pair_str]
            else:
                pred_col[i] = avg_pred
                print(pair_str+" NOT FOUND!")
        print("\n",end="")
        report_for_affinity(aff_col,pred_col)
    return


