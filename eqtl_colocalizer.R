#eqtl_colocalizer.R
#This script performs colocalization between a trait of the user's interest and the GTEx single trait eQLTs from
# the GTEx online GUI. It will perform colocalization for each gene-tissue pair in the GTEx csv file downloaded
# from the online GUI. NOTE it currently assumes you wish to use The GTEx_v7 data available on the voight lab LPC
# single trait eQLTs from at /project/voight_datasets/GTEx_v7/TissueSpecific.

#This script outputs to the directory it is run from and expects a "eQTL_config.R" file to be in the directory 

library(coloc)
library(data.table)

#Read in the arguments from the config file
source("eQTL_config.R")

#Print all of the config file settings to screen or the stnd out file
print("trait")
print(trait)
print("trait_file_header_info")
print(traitFilePath)
print(trait_A1col)
print(trait_A2col)
print(trait_SNPcol)
print(trait_CHRcol)
print(trait_BPcol)
print(trait_Pcol)
print(trait_Ncol)
print(trait_MAFcol)
print("traitType:")
print(traitType)
print(traitProp)
print(gtexCsvPath)
print("Locus Info:")
print(chrom)
print(colocStart)
print(colocStop)


#read in the tissue specific eQLT summary file with the file names added
trait_region = fread(file=traitFilePath, sep="\t", header=TRUE)
print("trait input file successfully loaded")

#add "trait" to the Allele fields, the MAF field, and the SNP field to avoid confusion when we compare to the eQTL alleles later 
trait_A1col_str = paste(trait_A1col,"_trait",sep="")
trait_A2col_str = paste(trait_A2col,"_trait",sep="")
trait_SNPcol_str = paste(trait_SNPcol,"_trait",sep="")
trait_MAFcol_str = paste(trait_MAFcol,"_trait",sep="")

colnames(trait_region)[colnames(trait_region)== trait_A1col] <- trait_A1col_str
colnames(trait_region)[colnames(trait_region)== trait_A2col] <- trait_A2col_str
colnames(trait_region)[colnames(trait_region)== trait_SNPcol] <- trait_SNPcol_str
colnames(trait_region)[colnames(trait_region)== trait_MAFcol] <- trait_MAFcol_str

#remove any alphabetical characters from the chromosome column
trait_region[[trait_CHRcol]] <- as.integer(gsub('[a-zA-Z]', '', trait_region[[trait_CHRcol]]))

str(trait_region)

#grab the snps that are within the start stop and on the correct chromosome from the trait file
trait_region = trait_region[trait_region[[trait_CHRcol]] == chrom & trait_region[[trait_BPcol]] >= colocStart & trait_region[[trait_BPcol]] <= colocStop,]

head(trait_region,3)

tissueTable = read.table(file="/project/voight_GWAS/wbone/eqtl_coloc/GTEx_Tissue_Summary_with_filenames.csv", sep=",", header=TRUE)
#Convert to strings from factors
tissueTable$Tissue = as.character(tissueTable$Tissue)
tissueTable$Filename = as.character(tissueTable$Filename)


#read in the csv file of eGene-Tissue pairs that came from the GTEx website
eGenes = read.table(file=gtexCsvPath, sep=",", header=TRUE)

eGenes$Gencode.Id = as.character(eGenes$Gencode.Id)
eGenes$Tissue = as.character(eGenes$Tissue)
eGenes$Gene.Symbol = as.character(eGenes$Gene.Symbol)

#loop through the eGene-Tissue pairs in eGenes and prep running COLOC
for(i in 1:nrow(eGenes)){

    geneID <- eGenes[i,"Gencode.Id"]
    gene <- eGenes[i,"Gene.Symbol"]
    tissue <- eGenes[i,"Tissue"]
    
    print(geneID)
    print(gene)
    print(tissue)

    # find the all pair file that contains the tissue of interest
    tissueLine <- tissueTable[tissueTable$Tissue == tissue,] 
    allpair_filename <- tissueLine$Filename

    #if the Filename field is NA then skip this eGene-Tissue pair
    if (is.na(allpair_filename)){
        print(tissue)
        print("This tissue is not available in the all pairs files currently")
        next
    }
    
    allpair_path = paste("/project/voight_datasets/GTEx_v7/TissueSpecific/", allpair_filename, sep="")
    eqtl_N <- tissueLine$NumberRNASeqandGTSamples

    #parentheses are causing issues too
    tissue_noSpace = gsub("\\(","",tissue)
    tissue_noSpace = gsub("\\)","",tissue_noSpace)
    # The tissue names have any whitespace in them and we want to use these in the output file names so replace " " with "_"
    tissue_noSpace = gsub("[[:space:]]","_",tissue_noSpace)
    
    # make a eGene-Tissue and trait prefix for file names
    out_prefix = paste(gene,geneID,tissue_noSpace,trait, sep="_")

    print("Grabing the all pairs data")    
    #generate the grep command for grabing the data from the chromosome of our locus
    grepCommand = paste0("| grep -P \"\t",chrom,"_\" | tr \"_\" \"\t\" > ",sep="")

    eGeneTissueInputFile = paste(gene,tissue_noSpace,chrom,colocStart,colocStop,".txt", sep="_" )
    # grab the eQTL region for this gene from the Tissue all pairs file
    #NOTE: We need to make the grep statement in a prior paste
    system(paste0("zcat ", allpair_path, grepCommand, eGeneTissueInputFile ))

    print("reading the all pairs data into R")
    #read the file we just generated from the grep command into R
    eGeneTissueInput = fread(file = eGeneTissueInputFile, sep="\t", header=FALSE)

    header = c("eGeneID","chr","bp","A1_eqtl","A2_eqtl","build","tss_dist", "ma_samples","ma_count","maf","pvalue_eQTL","slope","slope_se")

    colnames(eGeneTissueInput) = header

    print("Filtering the eGeneTissue_region file")
    #filter the input file to only include SNPs in the area we are interested in colocing 
    eGeneTissue_region = eGeneTissueInput[eGeneTissueInput$bp >= colocStart & eGeneTissueInput$bp <= colocStop,]

    print("Filtering on the geneID")
    #filter for the geneID of interest
    eGeneTissue_region = eGeneTissue_region[eGeneTissue_region$eGeneID == geneID,]

    if (nrow(eGeneTissue_region) == 0) {
      print("Warning: There wasn not an exact mathc on Ensembl ID. Likely this is due to a GTEX version  update.")
      #make a string that removes the everything after in the geneID
      noDecimalGeneID = gsub("\\..*","",geneID)
  
      #recreated eGeneTissue_region
      eGeneTissue_region = eGeneTissueInput[eGeneTissueInput$bp >= colocStart & eGeneTissueInput$bp <= colocStop,]
  
      #grep the simplified Ensembl ID
      possible_Ensembl_gene_lines <- eGeneTissue_region[grepl(noDecimalGeneID, eGeneTissue_region$eGeneID),]
  
      #check to make sure there is just one other Ensembl ID
      possible_Ensembl_genes <- unique(possible_Ensembl_gene_lines$eGeneID)
  
      if(length(possible_Ensembl_genes) == 1){
        print("Found a unique Ensembl ID so this analysis will continue using the Ensembl ID:")
        print(possible_Ensembl_genes)
    
        #this will be a single Ensembl ID string
        geneID = possible_Ensembl_genes
        eGeneTissue_region = eGeneTissue_region[eGeneTissue_region$eGeneID == geneID,]

      } else {
        print("The Ensembl ID from your GTEx csv was not able to be reliably mapped to an Enseml ID in the GTEx database, so this gene will be skipped:")
        print("geneID")
        next
      }    
    }

    print("adding rs numbers to the eQTL data")
    #add rs genegene,,numbers to the eGeneTissue_region DF
    #make a MarkerName column that is chr:pos
    eGeneTissue_region$MarkerName <- paste(eGeneTissue_region$chr,eGeneTissue_region$bp,sep=":")
  
    #merge the eGeneTissue_region and the dbSNP file with uniq IDs on MarkerName in order to add rs numbers
    #make uniq ID filename string 
    chromStr <- paste("chr",chrom,sep="")
    uniqID_filename = paste("/project/voight_GWAS/wbone/eqtl_coloc/",chromStr,"dbsnp_data_uniqID.txt",sep="")

    #read the uniqID filename into R
    uniqID_DF = fread(file = uniqID_filename , sep="\t", header=TRUE)

    print("merging the eqtl data with the dbSNP data on unique ID")
    #merge the files on the MarkerName
    eGeneTissue_region = merge(eGeneTissue_region, uniqID_DF, by = "MarkerName")

    #filter rows that dont have either A1_trait == A1_eqtl and A2_trait ==A2_eqtl OR A1_trait ==A2_eqtl and A2_trait == A1_eqtl
    eGeneTissue_region = eGeneTissue_region[eGeneTissue_region$ref == eGeneTissue_region$A1_eqtl & eGeneTissue_region$alt == eGeneTissue_region$A2_eqtl | eGeneTissue_region$alt == eGeneTissue_region$A1_eqtl & eGeneTissue_region$ref == eGeneTissue_region$A2_eqtl,]
    
    print("merging the trait and eqtl data on unique ID")
    #merge the trait and eGeneTissue region DFs on rs numbers
    colocInputFile = merge(eGeneTissue_region, trait_region, by.x="SNP", by.y=trait_SNPcol_str)

    #remove any NAs
    colocInputFile = colocInputFile[complete.cases(colocInputFile), ]

    #write colocInputFile to file for making locus zoom plots
    colocInputFile_outputStr = paste(out_prefix,"coloc_input_data.txt",sep="_")
    write.table(colocInputFile, file= colocInputFile_outputStr, sep="\t", row.names=FALSE, quote=FALSE)

    #remove any NAs
    #colocInputFile = colocInputFile[complete.cases(colocInputFile), ]

    print("Running coloc")
    #run coloc
    if (traitType == "cc"){

        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType, s=traitProp), dataset2=list(pvalues=colocInputFile$pvalue_eQTL, N=eqtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]])

    } else {

        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType), dataset2=list(pvalues=colocInputFile$pvalue_eQTL, N=eqtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]])

    }
    #prepare useful outputs
    coloc_results_summary = coloc_results$summary
    coloc_results_full = coloc_results$results

    #calculate pp4 / pp3 + pp4
    PP3andPP4 = coloc_results_summary[5] + coloc_results_summary[6]

    pp4_conditional = coloc_results_summary[6] / PP3andPP4

    #prep coloc output strings
    coloc_results_summary_outputStr = paste(out_prefix,"coloc_results_summary.txt",sep="_")
    coloc_results_full_outputStr = paste(out_prefix,"coloc_results_full.txt",sep="_")
    coloc_results_pp4_cond_outputStr = paste(out_prefix,"coloc_results_pp4_cond.txt",sep="_")

    #write to file
    write.table(coloc_results_summary, file=coloc_results_summary_outputStr, sep="\t", row.names=TRUE, quote=FALSE)
    write.table(coloc_results_full, file=coloc_results_full_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
    write.table(pp4_conditional, file=coloc_results_pp4_cond_outputStr, sep="\t", row.names=FALSE, quote=FALSE)




}


