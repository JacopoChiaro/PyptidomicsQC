import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import sys
import os
import upsetplot
from upsetplot import UpSet, plot
from mhcflurry import Class1PresentationPredictor
import io
import base64
#plt.style.use('ggplot')


## import replicates datasets ##
def PyptidomicsQC(list_of_dfs, HLAs):

    def create_file_dict_subset(list_of_dfs, minim=6, maxim=25):
        file_dict2 = {}
        new_list_of_dfs=[]
        for df in list_of_dfs:
            df = df.copy() 
            df2 = df[(df['Length'] > minim) & (df['Length'] <= maxim)]

            new_list_of_dfs.append(df2)
        
        return new_list_of_dfs

    def create_file_dict_9mers(list_of_dfs): ## modify to use list_of_dfs
        new_list_of_dfs = []
        for df in list_of_dfs:
            df2 = df[(df['Length'] == 9)]

            new_list_of_dfs.append(df2)
        
        return new_list_of_dfs

    def create_dict_lenght(x, replicate_df_Lenght): ## modify to use list_of_dfs
        dictionary = {}
        len_list = list(replicate_df_Lenght)
        
        for n in x:
            count=0
            for target in len_list:
                if n == target:
                    count+=1
            key = n
            value = count
            dictionary[key] = value

        return(dictionary)

    
    file_df_sub = create_file_dict_subset(list_of_dfs) ## future: to add possibility to change interval
    file_df_9mers = create_file_dict_9mers(list_of_dfs)


    ## Plot Peptide Lenght distribution ##
    def lenght_distr_plot(list_of_dfs, _title, **kwargs):
        col_names=[]
        for i,_ in enumerate(list_of_dfs, start=1):
            col_names.append(f'Replicate_{str(i)}')

        all_lenghts = []
        for df in list_of_dfs:
            all_lenghts = all_lenghts + list(df.Length)
        
        all_lenghts = list(set(all_lenghts))
        all_lenghts.sort()
        

        list_of_dict = []
        for df in list_of_dfs:
            dict_x = create_dict_lenght(all_lenghts, df.Length)
            list_of_dict.append(dict_x)

        list_of_lengths = []
        for i,d in enumerate(list_of_dict):
            value = [x for x in list(list_of_dict[i].values())]

            list_of_lengths.append(tuple(value))
        
        a = np.array(list_of_lengths)
        print(a.shape)
        arr = a.T
        print(arr.shape)

        list_of_perc = []
        for i,d in enumerate(list_of_dict):
            tot = sum(list(list_of_dict[i].values()))
            percent = [(x*100)/tot for x in list(list_of_dict[i].values())]
            
            list_of_perc.append(tuple(percent))
        
        a2 = np.array(list_of_perc)
        arr2 = a2.T

        hist_df = pd.DataFrame(data=arr, columns=col_names, index=all_lenghts)
        hist_percent_df = pd.DataFrame(data=arr2, columns=col_names, index=all_lenghts)
        
        return hist_df, hist_percent_df
    
    # #lenght_distr_plot(file_dict, f'{sample_name}_all_lenghts')
    #graph1, graph2 = lenght_distr_plot(file_df_sub, '8-25mers')
    hist_df, hist_percent_df = lenght_distr_plot(file_df_sub, '8-25mers')

    
    ## clean dataset from PTM ##
    def remove_PTM(pept):
        if '(' in pept:
            pept2 = re.sub(r"\(\+\d+\.\d+\)", "", pept)
            return pept2
        else:
            return pept

    def mass_spec_pept_cleanup(peptide):
        a = re.sub("^[A-Z]\.", "", peptide)
        if '.' in a:
            b = re.sub("\.[A-Z]$", "", a)
            return b
        else:
            return a

    ## https://upsetplot.readthedocs.io/en/latest/auto_examples/plot_missingness.html

    def Upset_plot(list_of_dfs, _title):    ## Fix the Upset plot code

        ## make a set of the total peptides found in the replicates

        merged_df = pd.concat(list_of_dfs)
        merged_df['No_PTM_peptides'] = [remove_PTM(x) for x in merged_df.Peptide]
        set_peptides = set(merged_df['No_PTM_peptides'])
        
        ## iterate through the replicates to make individual sets of peptides
        lists_of_peptides,replicates = [], []
        for x,df in enumerate(list_of_dfs, start=1):
            df['No_PTM_peptides'] = [remove_PTM(x) for x in df.Peptide]
            peptides = list(set(df['No_PTM_peptides']))
            lists_of_peptides.append(peptides)
            replicates.append(f'replicate_{str(x)}')

        df = pd.DataFrame(data=set_peptides, columns=['No_PTM_peptides'])

        for rep,plist in zip(replicates,lists_of_peptides):
            df[str(rep)] = [x in plist for x in df['No_PTM_peptides']]
        
        df['Counts'] = df[df.columns[1:]].sum(axis=1)
        df['Percent_Observed'] = [round(100*(x/len(replicates)),2) for x in df['Counts']]
        
        overlap_df = df.sort_values(by='Percent_Observed', ascending=False)
        
        upset_data = overlap_df.copy().drop(['No_PTM_peptides','Percent_Observed'], axis=1)

        upset_data2 = upset_data.set_index(list(upset_data.columns[:-1]))

        return overlap_df, upset_data2

    
    if len(list_of_dfs) > 1:
        overlap_df, upset_data2 = Upset_plot(list_of_dfs, f"{'8-25mers'}\n\n")
        ## modify the following
        #     upset_data2 = Upset_plot(file_dict, sample_name) ## modify these
        #     upset_data2 = Upset_plot(file_dict_sub, str(f'{sample_name}_8-20mers'))
        #     upset_data2 = Upset_plot(file_dict_9mers, str(f'{sample_name}_9mers'))
    else:
        overlap_df, upset_data2 = None, None

        




    ## MHC binding affinity prediction - MHC specificity Deconvolution##
    
    predictor = Class1PresentationPredictor.load()

    def create_counts_array(Netdf, hla_list):
        result=[]

        Netdf_sub = Netdf[Netdf['Length'] == 9]

        for elem in hla_list:
            Netdf_sub_sub = Netdf_sub[Netdf_sub['best_allele'] == elem]
            nB,WB,SB = 0,0,0
            ranks = list(Netdf_sub_sub['affinity_percentile'])
            for rank in ranks:
                if rank < 0.5:
                    SB+=1
                elif rank >2:
                    nB+=1
                else:
                    WB+=1
            result.append([nB,WB,SB])

        return np.array(result)



    def run_MHCaffinity_pred(list_of_dfs, allele_list):

        for df in list_of_dfs:
            df['No_PTM_peptides'] = [remove_PTM(x) for x in df.Peptide]

        MHC_runs=[]
        for df in list_of_dfs:
            df2 = df[df['Length'] < 16]
            run = predictor.predict(peptides=list(df2['No_PTM_peptides']), 
                                    alleles=allele_list, 
                                    include_affinity_percentile=True)
            run['Length'] = [len(x) for x in run.peptide]
            MHC_runs.append(run)

        count_arrays=[]
        percent_dfs=[]
        for run in MHC_runs:
            arr = create_counts_array(run, allele_list)
            res_df = pd.DataFrame(data=arr, columns=['nB','WB','SB'], index=[x for x in allele_list])
            count_arrays.append(res_df)


            tot = res_df['nB'].sum() + res_df['SB'].sum() + res_df['WB'].sum()

            res_df_percent = pd.DataFrame()
            res_df_percent['nB'] = [(x*100)/tot for x in res_df['nB']]
            res_df_percent['WB'] = [(x*100)/tot for x in res_df['WB']]
            res_df_percent['SB'] = [(x*100)/tot for x in res_df['SB']]
            res_df_percent.index = list(res_df.index)

            percent_dfs.append(res_df_percent)

        return MHC_runs, count_arrays, percent_dfs


    MHC_runs, count_arrays, percent_dfs = run_MHCaffinity_pred(list_of_dfs, HLAs)



    
    return hist_df, hist_percent_df, overlap_df, upset_data2, MHC_runs, count_arrays, percent_dfs