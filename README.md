# Test 2


#### Make a directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/
```
cd /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/
mkdir zsmith10
```
---
## FastQC v0.11.9 [Directory: 1_fastqc] : Assess Quality of Sequence Reads.
---
#### 1. Make a new `fastqc` directory and change into it.
```
mkdir 1_fastqc
cd 1_fastqc
```

#### 2. Symbolically link Solenopsis invicta sequence files to 1_fastqc directory.
```
ln -s ../../../raw_data/solenopsis_invicta_test2/*fastq .
```

#### 3. Load `fastqc` in `spack`.
```
spack load fastqc
```

#### 4. Process .fastq files into .fastqc files.
```
for file in *.fastq; do fastqc $file; done
```

#### 5. Secure copy files to your computer to view .html versions of .fastqc files.
##### 5a. Open a new local terminal.
```
cd desktop
mkdir test2_1_fastqc
cd test2_1_fastqc
scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/1_fastqc/*html' .
```
##### 5b. View files in browser to check quality, then return to the sphinx terminal.

---
## Skewer v0.2.2 [Directory: 2_skewer]: Trim Sequence Reads of Adapter Sequences.
---
#### 1. Make a new `skewer` directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10 and change into it.
```
cd ../
mkdir 2_skewer
cd 2_skewer
```

#### 2. Symbolically link .fastq files in the 2_skewer file.
```
ln -s ../../../raw_data/solenopsis_invicta_test2/*fastq .
```

#### 3. Search for adapter sequence counts in each file using the grep command. The -c flag will print counts of occurences.
##### Note: Illumina adapter sequence: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
```
grep -c 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG' *.fastq
```
###### Adapters found in each respective file for this test: 3459, 1106, 3105, and 964.

#### 4. Trim all .fastq files of synthetic adapter sequences using a for loop.
```
for f in *fastq; do BASE=$( basename $f | sed 's/_1.fastq//g'); echo $BASE;  /sphinx_local/software/skewer/skewer -t 2 -l 95 -Q 30 -x AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG $f -o $BASE; done
```
* the BASE=$( basename $f | sed 's/_1.fastq//g') function will define the basename of a file and remove the _1.fastq from all file names. Skewer will automatically add -trimmed.fastq to the base file names after analysis.
* The -t flag defines the number of threads for Skewer to use
* The -l flag defines minimum sequence lengths to keep
* The -Q flag defines the minimum Phred quality score for seqeuences

#### 5. Run `fastqc` on trimmed .fastq files to generate .fastqc files in html format.
```
fastqc *trimmed.fastq
```

#### 6. Secure copy trimmed .fastqc.html files to your home computer's terminal.
##### 6a. Open a new local terminal.
```
mkdir test2_2_skewer
cd test2_2_skewer
scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/2_skewer/*html' .
```
##### 6b. View files in browser to check quality and ensure sequences were removed, then return to the sphinx terminal.

---
## BWA v0.7.17 [Directory: 3_bwa]: Burrows-Wheeler Alignment (Aligning Sequences to a Reference)
---
#### 1. Create a directory for `bwa` analyses and change directories to it.
```
mkdir 3_bwa
cd 3_bwa
```

#### 2. Create a directory to house the reference genome and all pre-indexed files, then symbolically link them.
```
mkdir solenopsis_genome_index
cd solenopsis_genome_index
ln -s ../../../../raw_data/solenopsis_invicta/genome/UNIL_Sinv_3.0.fasta* .
```

#### 3. Return to the 3_bwa directory and symbolically link trimmed .fastq files from the 2_skewer directory.
```
cd ../
ln -s ../2_skewer/*-trimmed.fastq .
```

#### 4. Load `bwa` v0.7.17 and `samtools` v1.9 using `spack`
```
spack load bwa
spack load samtools@1.9%gcc@8.4.1
```

#### 5. Create a new script file (file extension = .sh) using `nano`. You may need to load `nano` using `spack`.
```
spack load nano
nano bwa.sh
```

#### 6. Here, we will use a for loop to: 1) define basenames for our .fastq files by removing the ".fastq" text and printing to the terminal using echo; 2) loading `bwa` and `samtools` using `spack`; and 3) running `bwa` to align S. invicta sequences to the reference genome, producing SAM files; and 4) piping bwa SAM file output into samtools to convert files into BAM format.
#### In the open nano file, enter the following script:
```
for file in *.fastq
do
    basename=$(echo "$file" | sed 's/.fastq//')
    echo $file
    echo $basename

    spack load bwa
    spack load samtools@1.9%gcc@8.4.1

    bwa mem -t 6 \
    solenopsis_genome_index/UNIL_Sinv_3.0.fasta \
    ${basename}.fastq \
    | samtools view -bSh \
    | samtools sort \
    -@ 3 -m 4G \
    -o ${basename}_sorted.bam
done
```
* the -t flag specified the number of threads for `bwa` to use.
* the -bSh flag in `samtools` converts large SAM files into smaller, binary BAM files with headers.
* the -@ flag designates threads for the job in `samtools`.
* the -m flag designates memory for the job in `samtools`.

#### 7. Run the script, which will produce sorted BAM files, using `bash` with the following command:
```
bash bwa.sh
```

#### 8. To assess number of mapped and supplemental reads, create a script for the following for loop to create read-mapping stats files with `samtools`. Use `bc`'s calculator tool to divide by 4 to correct the number of mapped reads
```
nano samtools_flagstat.sh
```
```
spack load samtools@1.9
spack load bc%gcc@8.4.1

for file in *_sorted.bam
do
        basename=$(echo $file | sed 's/-trimmed_sorted.bam//')
        samtools flagstat $basename-trimmed_sorted.bam > $basename-trimmed_sorted.stats
	
	echo $basename
        echo $(cat $basename-trimmed.fastq | wc -l)/4 | bc

done
```

#### 9. Add read groups to all BAM files for processing in GATK using `java` and `samtools`. First, load `java` (for Picard Tools) and `samtools`. 
```
spack load openjdk@11.0.8_10%gcc@8.4.1
spack load samtools@1.9
```

#### 10. Now, add read groups to each sample and index them with a for loop using Picard Tools.
```
nano picardtools.sh

for f in *_sorted.bam
do
        BASE=$( basename $f | sed 's/_sorted.bam*//g' )
        echo "BASE $BASE"
        
	java -jar /pickett_shared/software/picard-2.27.4/picard.jar \
		AddOrReplaceReadGroups \
		I=${BASE}_sorted.bam \
		O=${BASE}_sorted.RG.bam \
		RGSM=$BASE \
		RGLB=$BASE \
		RGPL=illumina \
		RGPU=$BASE
	samtools index ${BASE}_sorted.RG.bam
done
```  

#### 11. Exit the `nano` text file and run the script using `bash`.
```
bash picardtools.sh
```

#### 12. Double check that BAM files with read groups (.RG.bam files) actually had read groups added by checking the size of both .bam and their corresponding .RG.bam files. Add the flag `-lh` to view files in a list with human-readable sizes.
```
ls -lh *bam
```
#### 13. Correct read mapping percentages with `bc` tools by dividing by 4.
```
spack load bc%gcc@8.4.1
echo $(cat *-trimmed.fastq | wc -l)/4 | bc
```

---
## 4. GATK v4.1.8.1: Variant Calling [Directory: 4_gatk]
---
#### 1. First, return to your personal directory (/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10), make a new directory for your `gatk` analyses, and change into it.
```
cd ../
mkdir 4_gatk
cd 4_gatk
```
#### 2. Next, symbolically link the Solenopsis invicta reference genome, bam files with read groups (.RG.bam files) and index bam files with read groups (.RG.bam.bai files) to the 4_gatk directory.
```
ln -s /pickett_shared/teaching/EPP622_Fall2022/raw_data/solenopsis_invicta/genome/UNIL_Sinv_3.0.fasta .
ln -s ../3_bwa/*RG.bam* .
```
#### 3. Create a reference genome sequence dictionary (a .dict file) using `gatk` Note: `gatk` is natively installed on the server and is not loaded using `spack`.
```
/pickett_shared/software/gatk-4.2.6.1/gatk CreateSequenceDictionary -R UNIL_Sinv_3.0.fasta
```
#### 4. Index the reference genome .dict file using `samtools`, which create an indexed .fasta file with the .fasta.fai file type.
```
samtools faidx UNIL_Sinv_3.0.fasta
```

#### 5. Now, write a for loop to call SNPs and indels using the `gatk` HaplotypeCaller to compare variants. Files will be output in VCF and GVCF format. Standard VCF files only report detected variant calls, while GVCF files will allow us to recall SNPs/indels across samples by reporting all loci regardless of variants.
```
nano gatk_haplotypecaller.sh

for f in *_sorted.RG.bam
do
        BASE=$( basename $f | sed 's/-trimmed_sorted.RG.bam*//g' )
        echo "BASE $BASE"

        /pickett_shared/software/gatk-4.2.6.1/gatk HaplotypeCaller \
        -R UNIL_Sinv_3.0.fasta \
        -I ${BASE}-trimmed_sorted.RG.bam \
        -O ${BASE}.vcf

        /pickett_shared/software/gatk-4.2.6.1/gatk HaplotypeCaller \
        -R UNIL_Sinv_3.0.fasta \
        -I $f \
        -O ${BASE}.g.vcf \
        -ERC GVCF \
        -bamout ${BASE}_sorted.BG.realigned.bam
done

```

#### 6. Exit the `nano` text file and run the script using `bash`.
```
bash gatk_haplotypecaller.sh
```

#### 7. Copy the reference genome, reference genome index, all RG.bam and RG.bam.bai files, and all .vcf files to your computer using your local terminal.
```
scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/4_gatk/SRR*' .
scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/4_gatk/UNIL*' .
```

#### 8. View these files on IGV (igv.org). IGV allows you to visualize whether SNPs in the VCF files are supported, as well as visualize multiple read alignments.

---
## 5. GATK v4.1.8.1: Combining GVCFs and calling SNPs across sample groups [Directory: 5_gatk_combine]
---
#### 1. Make a new directory and change into it.
```
mkdir 5_gatk_combine
cd 5_gatk_combine
```

#### 2. Symbolically link all reference genome files.
```
ln -s /pickett_shared/teaching/EPP622_Fall2022/raw_data/solenopsis_invicta/genome/UNIL_Sinv_3.0.fasta
ln -s ../4_gatk/UNIL_Sinv_3.0.fasta.fai
ln -s ../4_gatk/UNIL_Sinv_3.0.dict .
```

#### 3. Link all VCF and GVCF files from the 4_gatk directory.
```
ln -s /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/4_gatk/*vcf .
```

#### 4. Write a script (`gatk_combine_gvcf_create.text.sh`) to generate a text command to combine GVCF files using the `gatk` CombineGVCF tool. 
*Define your reference as the Solenopsis invictus reference genome using the -R flag in gatk Combine GVCF tools. 
*Append each echo command to a new script file (gatk_combine_gvcf_run.sh) using `>>`. This will allow you to save the text output and rapidly regenerate the combined GVCF file. (Warning: extra spaces and tabs will break the script!)

```
nano gatk_combine_gvcf_create.text.sh
```
```
echo "/pickett_shared/software/gatk-4.2.6.1/gatk CombineGVCFs \\
-R UNIL_Sinv_3.0.fasta \\" >> gatk_combine_gvcf_run.sh

for f in *g.vcf
do
        echo "-variant $f \\" >> gatk_combine_gvcf_run.sh
done

echo "-O solenopsis_combined.g.vcf.gz" >> gatk_combine_gvcf_run.sh

```

#### 5. Now that you have created your `gatk_combined_gvcf_create.txt.sh` script, run it to create a script to generate a combined GVCF file.
```
bash gatk_combine_gvcf_create.text.sh
```

#### 6. This will create the aforementioned `gatk_combine_gvcf_run.sh` script. Run it to generate a combined GVCF file.
```
bash gatk_combine_gvcf_run.sh
```
#### 7. Now call variants again using the GenotypeGVCFs tool in `gatk` to generate a VCF file.
```
/pickett_shared/software/gatk-4.2.6.1/gatk --java-options "-Xmx10g" GenotypeGVCFs \
   -R UNIL_Sinv_3.0.fasta  \
   -V solenopsis_combined.g.vcf.gz \
   -O solenopsis_combined.vcf.gz
```

#### 8. Check the number of SNPs using `bcftools` to create a `.txt` file.
```
spack load bcftools
bcftools stats solenopsis_combined.vcf.gz > solenopsis_combined.vcf.stats.txt
```
*Variant information in file output appears so:
```
SN      0       number of samples:      4
SN      0       number of records:      59869
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 52557
SN      0       number of MNPs: 0
SN      0       number of indels:       7442
SN      0       number of others:       0
SN      0       number of multiallelic sites:   1235
SN      0       number of multiallelic SNP sites:       137
```

#### 9. Before filtering, unzip the `solenopsis_combined.vcf.gz` file in order to view files and ensure filters are working properly along the way.
```
gunzip solenopsis_combined.vcf.gz
```

#### 10. Filter for sequence depth and quality using manual parameters in the `gatk` VariantFiltration tool.
```
/pickett_shared/software/gatk-4.2.6.1/gatk VariantFiltration \
        -R UNIL_Sinv_3.0.fasta  \
        -V solenopsis_combined.vcf \
        -O solenopsis_combined.Basicfilters.vcf \
        -filter-name "basic_filter" -filter "QUAL > 20.0 && DP > 5"
```

#### 11. Check how many SNPs there are in the `Basicfilters` file using `bcftools` to create a `.txt` file. 
```
spack load bcftools
bcftools stats solenopsis_combined.Basicfilters.vcf > solenopsis_combined.Basicfilters.vcf.stats.txt
```

#### 12. Filter using the 'gatk' recommended filtering parameters to compare. `gatk` filters help with strand bias, filtering for quality based upon depth, etc. Note: There is no need to check SNP counts using `bcftools` because no information is actually removed by filtering, only flagged differently. Thus, output parameters, such as SNP counts, will be identical.
```
/pickett_shared/software/gatk-4.2.6.1/gatk VariantFiltration \
        -R UNIL_Sinv_3.0.fasta  \
        -V solenopsis_combined.vcf \
        -O solenopsis_combined.GATKfilters.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum" -filter "ReadPosRankSum < -8.0" 
```
#### 13. Gzip the following GVCF files:
*solenopsis_combined.vcf 
*solenopsis_combined.Basicfilters.vcf 
*solenopsis_combined.GATKfilters.vcf 
```
gzip solenopsis_combined.vcf
gzip solenopsis_combined.Basicfilters.vcf
gzip solenopsis_combined.GATKfilters.vcf
```

#### 13. Compare filtered and unfiltered files using `bcftools`.
```
bcftools stats -f "PASS,." solenopsis_combined.vcf.gz >solenopsis_combined.vcf.stats.txt
bcftools stats -f "PASS,." solenopsis_combined.Basicfilters.vcf.gz >solenopsis_combined.Basicfilters.vcf.stats.txt
bcftools stats -f "PASS,." solenopsis_combined.GATKfilters.vcf.gz >solenopsis_combined.GATKfilters.vcf.stats.txt
```

#### 13. Search for the number of SNPs in each `.stats.txt` file.
grep 'number of SNPs:' *stats.txt

SNP Output:
```
solenopsis_combined.Basicfilters.vcf.stats.txt:SN       0       number of SNPs: 3194
solenopsis_combined.GATKfilters.vcf.stats.txt:SN        0       number of SNPs: 36851
solenopsis_combined.vcf.stats.txt:SN    0       number of SNPs: 52557
```
Indel Output:
```
solenopsis_combined.Basicfilters.vcf.stats.txt:SN       0       number of indels:       333
solenopsis_combined.GATKfilters.vcf.stats.txt:SN        0       number of indels:       5130
solenopsis_combined.vcf.stats.txt:SN    0       number of indels:       7442
```
