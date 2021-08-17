#!/usr/bin/env bash
# Author Jonathan Lombardino, Burton Lab, University of Wisconsin Madison
# Convert NCBI GFF into BED file with strand orientation and base offset included
# Usage: get_intergenic.bash myfile.gff

gff=$1
name=${gff%.gff}
sorted_gff=""$name".sorted.gff"
gfile=""$name".genome"
sorted_gfile=""$name".sorted.genome"
bed_gff=""$name".bed" 
sorted_bed=""$name".sorted.bed"
intergenic=""$name"_intergenic.bed"
intergenic_stranded=""$name"_intergenic.stranded.bed"
all_bed=""$name"_all_features.bed"
all_bed_stranded=""$name"_all_features.stranded.bed"
all_gff=""$name"_all_features.gff"
all_gff_stranded=""$name"_all_features.stranded.gff"

## BED: 1:chrom,2:start,3:end,4:name,5:score,6:strand, ADDED: 7:source, 8:region
## GFF: 1:chrom,2:source,3:type,4:start,5:end,6:?,7:strand,8:offset,9:ID

echo "-------------| Formating and Sorting input GFF file"
### Sorting and Formatting input GFF file ###
awk '($3!="region") {print}' $gff | bedtools sort -i stdin > $sorted_gff

echo "-------------| Converting input GFF to BED formats and genome files"
### Create tab-delimited BED file skipping the genome entries (region): chrom, start(0),end(0-non-inclusive), name, score, strand ###
awk -F '\t' '{printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,int($4)-1,$5,$9,$8,$7,$2,$3);}' $sorted_gff > $bed_gff
awk -F '\t' '($3=="region") {printf("%s\t%d\n",$1,int($5)-1);}' $gff > $gfile  #make genome file using region entry that encompasses assembly.

echo "-------------| Getting intergenic regions from sorted GFF"
### Get the complement by strand and merge into one file at the end, code obtained from https://www.biostars.org/p/117925/  ###
awk '($7=="-") {print}' $sorted_gff | bedtools complement -i stdin -g $gfile | sed 's/$/\t-/' > $intergenic #get minus strand complement 
awk '($7=="+") {print}' $sorted_gff | bedtools complement -i stdin -g $gfile | sed 's/$/\t+/' >> $intergenic #get plus strand complement 

### Remove overlapping features by strand to get unique intergenic regions that do not overlap with coding sequences ###
awk -F '\t' '{printf("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n", $1,$2,$3,".",0,$4,".",".");}' $intergenic | \
bedtools sort -i stdin | bedtools subtract -s -a stdin -b $bed_gff | \
awk -F '\t' '{printf("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n", $1,$2,$3,"ID=inter"NR,0,$6,"RefSeq","intergenic");}' | \
cat - $bed_gff >  $all_bed_stranded #combine intergenic regions with original file

### Perform Bedtools complement but without any stranded information ###
bedtools complement -i $sorted_bed -g $sorted_gfile | 
awk -F '\t' '{printf("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n", $1,$2,$3,"ID=inter"NR, 0," ","Refseq","intergenic");}' | \
cat - $bed_gff | bedtools sort -i stdin > $all_bed   

echo "-------------| Converting BED outputs to GFF format"
### Convert combined bed back into GFF format ###
awk -F '\t' '{printf("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",$1,$7,$8,int($2)+1,$3,".",$6,$5,$4);}' $all_bed_stranded | bedtools sort -i stdin > $all_gff_stranded
awk -F '\t' '{printf("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",$1,$7,$8,int($2)+1,$3,".",$6,$5,$4);}' $all_bed | bedtools sort -i stdin > $all_gff


