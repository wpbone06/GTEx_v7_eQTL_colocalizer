#!/usr/bin/python2.7
#dbsnp_uniqID_maker.py

# uses some pandas commands to take in the output of grab_dbsnp150_chrs.sh and make a unique identifier. This is all done to make the rsID merger possible in eqtl_colocalizer.R
# designed to take the chr#_snp150_adbridged.txt files in as input and output all the unique ID files in the same directory


import pandas as pd


def main():
    chr_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]

    for chrom in chr_list:

        print chrom

        inputfileStr = chrom + "_snp150_adbridged.txt"

        #load the chromosome file into python as a pandas DF
        dbsnp_data = pd.read_table(inputfileStr,low_memory=False,sep="\t",header=None)
    
        #add a header
        dbsnp_data.columns = ["chrom", "chromStart", "chromEnd","SNP","score","strand","refNCBI","refUCSC","ref", "alt"]

        #strip the "chr" string from the chrom colum
        dbsnp_data['chrom'] = dbsnp_data['chrom'].map(lambda x: x.lstrip('chr'))

        #convert chromEnd to a string
        dbsnp_data['chromEnd'] = dbsnp_data['chromEnd'].apply(str)

        #make a MarkerName column
        dbsnp_data['MarkerName'] = dbsnp_data[['chrom', 'chromEnd']].apply(lambda x: ':'.join(x), axis=1)

        #write to file
        outfileStr = chrom + "dbsnp_data_uniqID.txt"
        dbsnp_data.to_csv(outfileStr,sep="\t",index=False)


if __name__ == "__main__":
    main()
