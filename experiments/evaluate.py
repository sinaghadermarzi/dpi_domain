#!/usr/bin/env python3

from sklearn.metrics import average_precision_score, roc_auc_score, mean_absolute_error,mean_squared_error, matthews_corrcoef,  f1_score, confusion_matrix, classification_report
from scipy.stats import spearmanr,pearsonr
from lifelines.utils import concordance_index
import numpy 

def concise(vals):
    if isinstance(vals, (int, float)):
        outstr = '{:.3f}'.format(vals)
    else:
        mean_val = numpy.mean(vals)
        std_val = numpy.std(vals)
        outstr = '{:.3f}'.format(mean_val)+"±"+'{:.3f}'.format(std_val)
    return outstr


def report_for_binary(tru_val,pred,print_report = True):
    pred_array= numpy.array(pred)
    tru_array = numpy.array(tru_val)
    auc = roc_auc_score(tru_array,pred_array)
    aupr = average_precision_score(tru_array,pred_array)
    rate_of_prediction= numpy.sum(tru_array)/len(tru_array)
    th = numpy.percentile(pred,(1-rate_of_prediction)*100)
    bin_pred = (pred_array>th).astype(int)
    cm = confusion_matrix(tru_array,bin_pred)
    mcc = matthews_corrcoef(tru_array,bin_pred)
    f1 = f1_score(tru_array,bin_pred)
    sensitivity = cm[0,0]/(cm[0,0]+cm[0,1])
    specificity = cm[1,1]/(cm[1,0]+cm[1,1])
    if print_report:
        print("AUC:\t",concise(auc))
        print("AUPR:\t",concise(aupr))
        print("MCC:\t",concise(mcc))
        print("F1:\t",concise(f1))
        print("SEN:\t",concise(sensitivity))
        print("SPEC:\t",concise(specificity))
    return auc,aupr,mcc,f1, sensitivity, specificity

def report_for_affinity(tru,pred,print_report = True):
    ci = concordance_index(tru,pred)
    mae = mean_absolute_error(tru,pred)
    spearman = spearmanr(tru,pred)[0]
    pearson = pearsonr(tru,pred)[0]
    mse = mean_squared_error(tru,pred)
    if print_report:
        print("CI:\t",concise(ci))
        print("MAE:\t",concise(mae))        
        print("MSE:\t",concise(mse))
        print("PEA:\t",concise(pearson))
        print("SP:\t",concise(spearman))
    return ci,mae,mse, pearson, spearman



if __name__== "__main__":
    report_for_affinity([3,4,5],[2.5,7.1,3.25])
    print("\n")
    report_for_binary([0,0,0,1,1],[0.1,0.2,0.5,0.3,0.4])
