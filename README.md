In this repository you will find all of the code to generate dependency files and run the eQTL colocalizer pipline. Given a GWAS signal of your chosing this pipeline will run COLOC on your signal and all the GTEx v7 single-tissue eQTLs.

Preparing the pipeline:

Download These files:
* The hg37 dbSNP flat file of your choosing
* The all SNPs tissue specific GTEx v7 files from the GTEx website: https://www.gtexportal.org/home/datasets

There are also some depency files you will need to generate:

* Make a dependency file directory and run grab_dbsnp150_chrs.sh on your hg37 dbSNP flat file.
* Next run dbsnp_uniqID_maker.py in the same directory to generate your dbSNP_data_uniqID files.
* Make sure you to have the GTEx_Tissue_Summary_with_filenames.csv in the depency files directory, and to update the path to this file in eqtl_colocalizer.R


Instructions on running the eQTL colocalizer pipeline:

- Search your lead SNP on the GTEx website: https://www.gtexportal.org/home, and download the CSV file of the single-tissue-eQTLs.

- Create an analysis directory and add your eqtl_colocalizer.R, your GTEx CSV file and eqtl_config.R file corresponding the GWAS signal you would like to perform eQTL colocalization analysis on. Then execute the eqtl_colocalizer.R in that directory (eg Rscript ./eqtl_colocalizer.R)
