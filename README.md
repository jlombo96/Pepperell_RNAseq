# Pepperell_RNAseq
Code to perform intergenic RNA-seq analysis 
## Intergenic Analysis Pipeline 

#### Step 1: Prepare files

Requirements: bedtools v2.29.2 (though anything above v2 will *likely* work.) 

Run the code below to generate the following files
```
bash get_intergenic.bash myintputgenome.gff
```

The following files should be made
|File	|Description	|
|---	|---	|
|.sorted.gff	|	sorted input GFF	|
|.genome	|	bedtools genome file describing chromosome sizes|
|.bed	|	GFF file in BED format	|
|.intergenic.bed	|	intergenic regions from bedtools complement (non-stranded)|
|.intergenic.stranded.bed	|	intergenic regions from bedtools complement(stranded)|
|_all_features.bed	|	combined input GFF and intergenic regions in BED format|
|_all_features.stranded.bed	|	" " but stranded|
|_all_features.gff	|	combined input GFF and intergenic regions in GFF format|
|_all_features.stranded.gff	|	" " but stranded|

Try running this on your server to confirm bedtools is up-to-date and working. Otherwise I'll upload the files onto this repository. 

#### Step 2: Make Windows (optional):

Run the following code to generate windows from the output .sorted.genome from **Step 1**
```
bash makewindows.bash myinputgenome.sorted.bed windowsize windowstep
```
A file titled myinputgenome_windows_size{}_step{}.bed will be generated. This can be used to map RNA-seq reads to a sliding window in **Step 3**

#### Step 3: Map RNA-seq reads to intervals via bedtools (for visual, preliminary results)

Place the mapcoverage.bash in a folder containing your RNA-seq BAM files. **NOTE:** The following code may require the input BAM files to be sorted. You can do so with samtools sort.

```
bash mapcoverage.bash myintervalfile.{bed/gff/vcf}
```
Here you can use a feature file generated in **Step 1** (i.e. the .all_fatures.stranded.bed) or the sliding window file generated in **Step 2**. The script will loop through each of the sorted BAM files and report the read-depth across each interval in the provided file. The final column includes the name of the BAM file that generated the corresponding data (which can be pivoted to columns in R for downstream summary statistics). 

The current bedtools coverage parameters will output the following results for each interval:
1. The number of features in B that overlapped the A interval.
2. The number of bases in A that had non-zero coverage.
3. The length of the entry in A.
4. The fraction of bases in A that had non-zero coverage.

NOTE: the current code requires same-strandedness (-s) and treats split reads as a single interval (no -split flag). If you would like these changed let me know. See [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) for more info.




 

