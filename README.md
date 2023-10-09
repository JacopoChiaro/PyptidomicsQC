![PyptidomicsQC_logo](https://github.com/JacopoChiaro/PyptidomicsQC/assets/70270016/21b50798-8dad-40ff-9e11-d8130b372cb2)

# PyptidomicsQC

PyptidOmicsQC is a web application developed in Flask, offering a user-friendly interface, and facilitating quick visualization of QC graphs through a dashboard. The script was written in Python, specifically using version 3.9.7 (from Anaconda version 2021.11), and it operates entirely within a Python environment.
Briefly, this tool automatically provides i) peptide length distribution (in raw numbers and percentage) of all the provided replicates at once side by side. ii) Representation of overlap between replicates is performed using the UpSet plot (from Upsetplot package: https://upsetplot.readthedocs.io/en/stable/), run with default settings. iii) Stacked bar plot deconvoluting the eluted peptides' allele specificity. Peptide specificity was assessed by MHCflurry package using allele rank for the best allele predicted. Peptides ranking below 0.5% were considered “strong binders”, while peptides ranking between 0.5% and 2% were considered “weak binders”. Peptides scoring above 2% are considered “non-binders”.

# Setup and run PyptidomicsQC

1. Download repository
2. install the rquired packages listed in "requirements.txt"
3. run "app.py" using flask and copy the obtained link in a browser.
4. to test the software use the example files contained in the "input folder"
5. upload the input files following the instruction shown on the web page.
6. check that the obtained output matches witht the one present in the "output" folder.

In case of problems don't hesitate to contact us at jacopo.chiaro@helsinki.fi or vincenzo.cerullo@helsinki.fi
