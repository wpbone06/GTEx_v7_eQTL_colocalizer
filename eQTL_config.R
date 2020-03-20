#ID for user's trait of interest. (Can be any string)
trait = "AD"
#path to the input files
traitFilePath = "/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
#column IDs from trait file
trait_A1col = "A1"
trait_A2col = "A2"
trait_SNPcol = "SNP"
trait_CHRcol = "CHR"
trait_BPcol = "BP"
trait_Pcol = "P"
trait_Ncol = "Nsum"
trait_MAFcol = "MAF"
#trait info not in the input file
#traitType is set either to "cc" or "quant"
traitType = "cc"
#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
traitProp = 0.157888494
gtexCsvPath = "/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/rs4308_GTEx_COLOC_AD/rs4308_GTEx_Portal_eQTLs_2020-01-07.csv"
#locus information for running coloc. Currently this assumes these genomic positions to be from hg19
chrom = 17
colocStart = 61309625
colocStop = 61809625
 
