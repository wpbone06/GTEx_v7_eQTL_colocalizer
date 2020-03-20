#grab_dbsnp150_chrs.sh
# This script splits the /project/voight_datasets/dbSNP/snp150.txt into chromosomes, removes triallelic loci, and removes extra fields to allow for rs numbers to be added to the eGene files easily


#$1 is the list of chromosomes you wish to grep out of the dbsnp150.txt file


while read chr
do
    echo $chr

    #make the grep string we need to search
    chromGrep=$chr"\t"
    echo $chromGrep

    #make output file string
    chrom_dbsnp_file=$chr"_snp150_adbridged.txt"
    echo $chrom_dbsnp_file


    #perform grep
    grep -P $chromGrep /project/voight_datasets/dbSNP/snp150.txt | grep -P "genomic\tsingle" | grep -v TriAllelic | grep -v -P "/./" | cut -f 2-10 | tr "/" "\t" > ./$chrom_dbsnp_file

done < $1

