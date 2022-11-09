import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import sys
import os
import upsetplot
from upsetplot import plot
from mhcflurry import Class1PresentationPredictor
plt.style.use('ggplot')



## input section

sample_name = input('Type sample name: ')

input_str,paths,count = '.',[],1

while input_str != '':
    input_str = input(f'insert a path of replicate {count}: ')
    count+=1
    if input_str != '':
        paths.append(input_str)

alleles = input('write down allele list for MHC binding affinity prediction with comma (no spaces): ')

Alleles_list = alleles.split(',')
print(Alleles_list)

## create a folder to collect all the outputs
home = '/Users/jacopoch'
proj_path = (f'{home}/Desktop/{sample_name}')
os.makedirs(proj_path)



## import replicates datasets ##

def create_file_dictionary(list_of_filepaths):
    file_dict = {}
    for i,file in enumerate(list_of_filepaths):
        key = i
        df = pd.read_csv(file)

        file_dict[key] = df
    
    return file_dict

def create_file_dict_subset(file_dict, minim=6, maxim=20):
    file_dict2 = {}
    for i,_ in enumerate(file_dict):
        key = i
        df = file_dict[i].copy() 
        df = df[(df['Length'] > minim) & (df['Length'] <= maxim)]

        file_dict2[key] = df
    
    return file_dict2

def create_file_dict_9mers(file_dict):
    file_dict3 = {}
    for i,_ in enumerate(file_dict):
        key = i
        df = file_dict[i].copy() 
        df = df[(df['Length'] == 9)]

        file_dict3[key] = df
    
    return file_dict3

def create_dict_lenght(x3, replicate_df_Lenght):
    dictionary = {}
    len_list = list(replicate_df_Lenght)
    
    for n in x3:
        count=0
        for target in len_list:
            if n == target:
                count+=1
        key = n
        value = count
        dictionary[key] = value

    return(dictionary)

file_dict = create_file_dictionary(paths)
file_dict_sub = create_file_dict_subset(file_dict) ## future: to add possibility to change interval
file_dict_9mers = create_file_dict_9mers(file_dict)

## Plot Peptide Lenght distribution ##
def lenght_distr_plot2(file_dict, sample_name, **kwargs):
    
    all_lenghts = []
    for i,_ in enumerate(file_dict):
        all_lenghts = all_lenghts + list(file_dict[i].Length)
    
    all_lenghts = list(set(all_lenghts))
    all_lenghts.sort()
    
    list_of_dict = []
    for i,_ in enumerate(file_dict):
        dict_x = create_dict_lenght(all_lenghts, file_dict[i].Length)
        #print(dict_x)
        list_of_dict.append(dict_x)
    
    #print(list_of_dict)

    list_of_perc = []
    for i,d in enumerate(list_of_dict):
        tot = sum(list(list_of_dict[i].values()))
        percent = [(x*100)/tot for x in list(list_of_dict[i].values())]
        
        list_of_perc.append(tuple(percent))
    
    a = np.array(list_of_perc)
    arr = a.T
    #print(arr)
    
    col_names=[]
    for i,_ in enumerate(file_dict):
        col_names.append(f'Replicate_{str(i+1)}')
    #print(col_names)
        
    hist_percent_df = pd.DataFrame(data=arr, columns=col_names, index=all_lenghts)
    

    hist_percent_df.plot(kind='bar', edgecolor='black', alpha=0.7, rot=0)
    plt.xlabel("Peptide Lenght")
    plt.ylabel("Percentage of unique peptides")
    plt.title(sample_name)
    plt.tight_layout()
    plt.savefig(f'{proj_path}/{sample_name}_peptide_lenght_distribution.png', dpi=300, format='png')
    plt.show()

lenght_distr_plot2(file_dict, f'{sample_name}_all_lenghts')
lenght_distr_plot2(file_dict_sub, f'{sample_name}_8-20mers')

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

def Upset_plot(file_dict, sample_name):
    files=[]
    for x,_ in enumerate(file_dict):
        files.append(file_dict[x])    

    merged_df = pd.concat(files)
    merged_df['No_PTM_peptides'] = [remove_PTM(x) for x in merged_df.Peptide]
    
    lists_of_peptides,replicates = [], []

    for file in files:
        file['No_PTM_peptides'] = [remove_PTM(x) for x in file.Peptide]

    for x,file in enumerate(files):
        peptidies = list(set(file['No_PTM_peptides']))
        lists_of_peptides.append(peptidies)
        replicates.append(f'replicate_{str(x+1)}')

        
    set_peptides = set(merged_df['No_PTM_peptides'])

    df = pd.DataFrame(data=set_peptides, columns=['No_PTM_peptides'])


    for rep,plist in zip(replicates,lists_of_peptides):
        df[str(rep)] = [x in plist for x in df['No_PTM_peptides']]
    
    df['Counts'] = df[df.columns[1:]].sum(axis=1)
    df['Percent_Observed'] = [round(100*(x/len(replicates)),2) for x in df['Counts']]
    
    df = df.sort_values(by='Percent_Observed', ascending=False)
    
    upset_data = df.copy().drop(['No_PTM_peptides','Percent_Observed'], axis=1)

    upset_data2 = upset_data.set_index(list(upset_data.columns[:-1]))
    
    
    plot(upset_data2,element_size=None,show_counts=True, show_percentages=True)
    plt.title(sample_name)
    plt.savefig(f'{proj_path}/{sample_name}_replicates_overlap.png', dpi=300, format='png')
    plt.show()

if len(file_dict_sub) > 1:
    Upset_plot(file_dict, sample_name)
    Upset_plot(file_dict_sub, str(f'{sample_name}_8-20mers'))
    Upset_plot(file_dict_9mers, str(f'{sample_name}_9mers'))
else:
    pass



## MHC binding affinity prediction - MHC specificity Deconvolution##
print('binding affinity prediction start...\n\n\n')
predictor = Class1PresentationPredictor.load()

def create_counts_array(Netdf, hla_list):
    result=[]

    Netdf_sub = Netdf[Netdf['Lenght'] == 9]

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

def plot_clustered_stacked3(dfall, sample_name, labels=None, title="HLA-specificity deconvolution",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    plt.figure(figsize=(10,4))
    axe = plt.subplot(111)
    
    
    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      alpha=0.7,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col)-(1/(n_col*2)))
                
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width((1 / float(n_df+1)))

    axe.set_xticks(((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)-(1/((n_col*n_df)-0.5)))
    axe.set_xticklabels(df.index, rotation = 0)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))
    
    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)
    
    plt.ylabel("Percentage of specific peptides")
    plt.tight_layout()
    plt.savefig(f'{proj_path}/{sample_name}_{title}.png', dpi=300, format='png')
    plt.show()
    
    return axe
    

def run_MHCaffinity_pred(file_dict_9mers, allele_list, sample_name, proj_path):
    files=[]
    for x,_ in enumerate(file_dict_9mers):
        files.append(file_dict_9mers[x])    

    lists_of_peptides,replicates = [], []

    for file in files:
        file['No_PTM_peptides'] = [remove_PTM(x) for x in file.Peptide]
    
    MHC_runs=[]
    for x,file in enumerate(files):
        run = predictor.predict(peptides=list(file['No_PTM_peptides']), 
                                alleles=allele_list, 
                                include_affinity_percentile=True)
        run['Lenght'] = [len(x) for x in run.peptide]
        run.to_csv(f'{proj_path}/{sample_name}_MHC_binding_rep_{str(x+1)}.csv', index=False)
        MHC_runs.append(run)
        
    count_arrays=[]
    percet_dfs=[]
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
        
        percet_dfs.append(res_df_percent)
        
    return count_arrays, percet_dfs


count_arrays, percet_dfs = run_MHCaffinity_pred(file_dict_9mers, Alleles_list, sample_name, proj_path)


plot_clustered_stacked3(count_arrays, 
                        sample_name, labels=[f"Replicate {str(x+1)}" for x,_ in enumerate(count_arrays)], 
                        title="9mers HLA-specificity deconvolution")


plot_clustered_stacked3(percet_dfs, 
                        sample_name, labels=[f"Replicate {str(x+1)}" for x,_ in enumerate(percet_dfs)], 
                        title="9mers HLA-specificity deconvolution percent")