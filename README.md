# AMbeRland


## Data: _in silico_ minimum inhibitory concentrations (MIC) 

We obtained whole genome sequence data by searching the NCBI Assembly database for the search term "Klebsiella pneumoniae", utilising the rentrez R package [1], and for each resulting entry fetching the related GenBank sequence. Sampling location and date was found from linking each assembly entry to its related NCBI BioSample entry, discarding those without location data. We then employed  an XGBoost-based machine learning package [2, 3] to predicts MICs _in silico_ for each antibiotic, for each genome, in each country, for every year of data available. 

User can provide custom data file, which has to be in the following tabulated format: "ID”, “Country", "Year", "Antibiotic”, ”MIC”.


## Visualisation

* Population distribution by year;
* Principal Component Analysis (PCA) [5]. 
* t-Distributed Stochastic Neighbor Embedding (t-SNE)  analysis [6].

For population distribution plots, we have log2(MIC) values on x-axis and collection year on y-axis, with gray circles showing all antibiotic concentration-year combinations present in the dataset and the size of black circle identifying  a proportion of  isolates inhibited by the specific level of antimicrobial agent. The closer a radius of black circle is to radius of gray circle, the closer this proportion is to one. We utilize color as a third dimension to visualise the log2(MIC) values per each isolate-antibiotic combination for PCA & t-SNE visualisation.

## Extension

Comparing [AMR surveillance data](https://amr.theodi.org/programmes/atlas) and _in silico_ MICs for _Klebsiella pneumoniae_:
[rretkute.shinyapps.io/wellcome_data_reuse](https://rretkute.shinyapps.io/wellcome_data_reuse/)

## References
[1] Winter  (2017) rentrez: an R package for the NCBI eUtils API The R Journal 9(2):520-526
[2] Nguyen et al. (2018).   Developing an in silico minimum inhibitory concentration panel test for Klebsiella pneumoniae. Scientific Reports volume 8, Article number: 421, DOI: 10.1038/s41598-017-18972-w
[3] https://github.com/PATRIC3/mic_prediction
[4]   Suzuki et al. (2011) Antimicrobial-susceptible patterns of Staphylococcus aureus isolated from surgical infections: a new approach.   J Infect Chemother.  17(1):34-9. DOI: 10.1007/s10156-010-0096-y
[5]  van der Maaten  (2014)  Accelerating t-SNE using Tree-Based Algorithms.  Journal of Machine Learning Research 3221-324. 