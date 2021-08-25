# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 08:54:19 2021

@author: moran
"""


'''
1. read FC file
2. read process file
3. TiPA function
4. create matrix + GO term2 name function
'''

import csv
import numpy as np
import pandas as pd
from scipy import stats
#import os


tissues = ['Adipose - Subcutaneous',
 'Adipose - Visceral (Omentum)',
 'Artery - Aorta',
 'Artery - Coronary',
 'Artery - Tibial',
 'Brain0',
 'Brain1',
 'Brain2',
 'Breast - Mammary Tissue',
 'Colon - Sigmoid',
 'Esophagus - Gastroesophageal Junction',
 'Esophagus - Mucosa',
 'Esophagus - Muscularis',
 'Heart - Atrial Appendage',
 'Heart - Left Ventricle',
 'Liver',
 'Lung',
 'Muscle - Skeletal',
 'Nerve - Tibial',
 'Ovary',
 'Pituitary',
 'Prostate',
 'Skin - Not Sun Exposed (Suprapubic)',
 'Skin - Sun Exposed (Lowerleg)',
 'Testis',
 'Thyroid',
 'Uterus',
 'Vagina',
 'Whole Blood']



# Create process and genes dictionary
# =====================================

processes_and_genes_file = csv.reader(open('process_and_associated_genes.csv', 'r'))  # processes and associated genes after filtering by GO domain ('biological process') and size 
proc_genes_dict = {}
for row in processes_and_genes_file:
    proc_genes_dict[row[0]] = row[1:]
    

def FC_tissue_v7(current_tissue):
    genes_FC_in_tissue_v7 = {}
    tis_file = csv.reader(open('tissues_FC_directory/' + current_tissue + '.csv', 'r'))  
    next(tis_file)
    for row in tis_file:
        gene = row[1]
        val = float(row[2])
        genes_FC_in_tissue_v7[gene] = val

    return genes_FC_in_tissue_v7


def TiPA_for_tissue(proc_genes_dict, current_tissue, genes_FC_tissue):
    
    ''' This function gets dictionary of processes and associated genes, 
    query tissue and FC vlaues for the tissue, and returns TiPA scores 
    alnog the size of the process'''
        
    scoring_and_length_per_tissue = {}

    for process in proc_genes_dict:
        list_of_tis_values = []
        for gene in proc_genes_dict[process]:
            if gene in genes_FC_tissue.keys():
                gene_value = float(genes_FC_tissue[gene])
                list_of_tis_values.append(gene_value)
        if len(list_of_tis_values) > 0:
            sorted_values = sorted(list_of_tis_values)
            tipa_score = stats.trim_mean(sorted_values, 0.1)  # calculation of TiPA score
            scoring_and_length_per_tissue[process] = (tipa_score)

    return scoring_and_length_per_tissue


def tipa_for_all_tissues():

    scoring_all_tissues = {}
    for tissue in tissues:
        genes_FC_tissue = FC_tissue_v7(tissue)
        scoring_all_tissues[tissue] = TiPA_for_tissue(proc_genes_dict, tissue, genes_FC_tissue)
    
    tipa_scores = {}
    for tissue in scoring_all_tissues:
        for process in scoring_all_tissues[tissue]:
            if process not in tipa_scores:
                tipa_scores[process] = {}
            tipa_scores[process][tissue] = scoring_all_tissues[tissue][process]
            
    scores_matrix = pd.DataFrame.from_dict(tipa_scores).transpose()
    scores_matrix.to_csv('tipa_scores_matrix_test_test_test.csv')

    return tipa_scores  



tipa_for_all_tissues()



    
