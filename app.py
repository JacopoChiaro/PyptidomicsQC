from flask import Flask, request, render_template, send_file, url_for
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style as style
import io
import re
import os
import base64
import upsetplot
from datetime import datetime
import tempfile
import shutil
from pyptidomicsQC_adapt import PyptidomicsQC
style.use('tableau-colorblind10')

app = Flask(__name__)
app.secret_key='PyptidomicsQC'

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():

    # Load the input files
    sample_name = request.form['sample_name']

    hlas = request.form['text']
    hlas = [x.strip() for x in hlas.split(',')]

    dfs = []

    file1 = request.files['file1']
    df1 = pd.read_csv(file1)
    dfs.append(df1)
    
    file2 = request.files['file2']
    if file2:
        df2 = pd.read_csv(file2)
        dfs.append(df2)
    else:
        df2 = None
    
    file3 = request.files.get('file3')
    if file3:
        df3 = pd.read_csv(file3)
        dfs.append(df3)
    else:
        df3 = None
    
    file4 = request.files.get('file4')
    if file4:
        df4 = pd.read_csv(file4)
        dfs.append(df4)
    else:
        df4 = None
    
    file5 = request.files.get('file5')
    if file5:
        df5 = pd.read_csv(file5)
        dfs.append(df5)
    else:
        df5 = None

    file6 = request.files.get('file6')
    if file6:
        df6 = pd.read_csv(file6)
        dfs.append(df6)
    else:
        df6 = None
    
    hist_df, hist_percent_df, overlap_df, upset_data2, MHC_runs, count_arrays, percent_dfs = PyptidomicsQC(dfs, hlas)


    ## Plot images
    ## Figure 1 ## Peptide Lengths distribution - Number of Peptides
    fig1, ax = plt.subplots()
    hist_df.plot(kind='bar', edgecolor='black', alpha=0.7, rot=0, ax=ax)
    ax.set_xlabel("Peptide Length")
    ax.set_ylabel("Number of unique peptides")
    ax.set_title('8-25mers')
    graph1 = plot_to_html(fig1)

    ## Figure 2 ## Peptide Lengths distribution - Percentage of Peptides
    fig2, ax = plt.subplots()
    hist_percent_df.plot(kind='bar', edgecolor='black', alpha=0.7, rot=0, ax=ax)
    ax.set_xlabel("Peptide Length")
    ax.set_ylabel("Percentage of unique peptides")
    ax.set_title('8-25mers')
    graph2 = plot_to_html(fig2)

    ## Figure 3 ## UpSet plot - Peptide Overlap
    fig3, ax = plt.subplots(figsize=(10,6))
    upsetplot.plot(upset_data2, element_size=None, show_counts=True, show_percentages=True, fig=fig3)
    ax.set_title(f"{'8-25mers'}\n")
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    graph3 = plot_to_html(fig3)

    ## Figure 4 ## Peptide MHC Binding Affinity Deconvolution Stacked
    fig4 = plot_clustered_stacked(count_arrays, 
                            'Number of specific peptides', labels=[f"Replicate {str(x)}" for x,_ in enumerate(count_arrays, start=1)], 
                            title="9mers HLA-specificity deconvolution")
    graph4 = plot_to_html(fig4)

    ## Figure 5 ## Peptide MHC Binding Affinity Deconvolution Stacked
    fig5 = plot_clustered_stacked(percent_dfs, 
                            'Percentage of specific peptides', labels=[f"Replicate {str(x)}" for x,_ in enumerate(percent_dfs, start=1)], 
                            title="9mers HLA-specificity deconvolution percent")
    graph5 = plot_to_html(fig5)


    ## process data for saving

    ## process MHC binding affinity runs
    MHCruns_to_merge = []
    for df in MHC_runs:
        df2=df[['peptide','affinity','best_allele','affinity_percentile','Length']]
        MHCruns_to_merge.append(df2)

    concat_MHC_runs = pd.concat(MHCruns_to_merge)
    concat_MHC_runs2 = concat_MHC_runs.drop_duplicates()

    ## process MHC binding affinity counts
    count_arrays_to_merge=[]
    for i,array in enumerate(count_arrays, start=1):
        df = pd.DataFrame(data=array)
        df['Replicate'] = [f'Replicate {i}' for x in df['nB']]
        count_arrays_to_merge.append(df)
        
    merged_counts = pd.concat(count_arrays_to_merge)

    ## process MHC binding affinity percentages
    percent_dfs_to_merge=[]
    for i,array in enumerate(percent_dfs, start=1):
        df = pd.DataFrame(data=array)
        df['Replicate'] = [f'Replicate {i}' for x in df['nB']]
        percent_dfs_to_merge.append(df)
        
    merged_percent = pd.concat(percent_dfs_to_merge)
    
    # Create a folder to save all the log files
    dt = str(datetime.now()).split('.')[0]
    dt2 = re.sub(':','-',dt)
    dt3 = re.sub(' ','-',dt2)
    folder_name = f'{sample_name}_{str(dt3)}'

    directory = os.getcwd()
    path = os.path.join(directory, f'output/{folder_name}')
    os.mkdir(path)

    # upload files in the newly created folder:
    hist_df.to_csv(f'{path}/peptide_length_distribution.csv')
    hist_percent_df.to_csv(f'{path}/peptide_percent_length_distribution.csv')
    overlap_df.to_csv(f'{path}/overlap_df.csv')
    upset_data2.to_csv(f'{path}/data_for_upsetplot.csv')
    concat_MHC_runs2.to_csv(f'{path}/peptide_binding_affinity.csv')
    merged_counts.to_csv(f'{path}/binders_count_array.csv')
    merged_percent.to_csv(f'{path}/binders_percentage_array.csv')


    # Render the template with the graphs
    return render_template('analyze.html', 
                           text=hlas, 
                           graph1=graph1, graph2=graph2, graph3=graph3, graph4=graph4, graph5=graph5, 
                           
                           folder_name=folder_name)                    


def make_temp_dir_and_download(folder_to_save):
    with tempfile.TemporaryDirectory() as tempdir:

        # compress folder into tempdir
        print(os.listdir(tempdir))
        shutil.move(f'output/{folder_to_save}', tempdir)
        print(os.listdir(tempdir))
        
        path_to_folder = os.path.join(tempdir, folder_to_save)
        shutil.make_archive(path_to_folder, 'zip', root_dir=tempdir, base_dir=folder_to_save)
        
        print(f'This is the tmp folder: {tempdir}')
        print(os.listdir(tempdir))

        path=os.path.join(tempdir, f'{folder_to_save}.zip')
        
        return send_file(path, mimetype='zip', as_attachment=True)

@app.route('/download/<filename>', methods=['GET', 'POST'])
def download(filename):
    print(f'\n\n\n\n{filename}    \n\n\n\n')
    # response = make_temp_dir_and_download(filename)       

    # return response
    shutil.make_archive(f'output/{filename}', 'zip', root_dir='output', base_dir=filename)
    return send_file(f'output/{filename}.zip', mimetype='zip', as_attachment=True)





def plot_to_html(fig):
    # Convert a Matplotlib plot to HTML for rendering in Flask
    img = io.BytesIO()
    print(f'\n\n\n{type(fig)}\n\n\n')
    fig.savefig(img, format='png')
    img.seek(0)
    plot_data = base64.b64encode(img.getvalue()).decode()
    return '<img src="data:image/png;base64,{}">'.format(plot_data)

def plot_clustered_stacked(dfall, y_lable_text, labels=None, title="HLA-specificity deconvolution",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
    labels is a list of the names of the dataframe, used for the legend
    title is a string for the title of the plot
    H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    colors = ['#ABABAB','#FF800E','#006BA4']
    fig = plt.figure(figsize=(8,4))
    axe = plt.subplot(111)
    
    
    for df in dfall: # for each data frame
        axe = df.plot(kind="bar",
                    linewidth=1,
                    color=colors,
                    edgecolor='black',
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
    
    plt.ylabel(y_lable_text)
    plt.tight_layout()
    return fig



if __name__ == '__main__':
    app.run(debug=True)