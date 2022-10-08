# Test 2


#### Make a directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/
```
cd /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/
mkdir zsmith10
```
---
## FastQC v0.11.9 [Directory: 1_fastqc] : Assess Quality of Sequence Reads.
---
#### 1. Make a new fastqc directory.
```
mkdir 1_fastqc
cd 1_fastqc
```

#### 2. Symbolically link Solenopsis invicta sequence files to 1_fastqc directory.
```
ln -s ../../../raw_data/solenopsis_invicta_test2/*fastq .
```

#### 3. Load fastqc in Spack.
```
spack load fastqc
```

#### 4. Process .fastq files into .fastqc files.
```
for file in *.fastq; do fastqc $file; done
```

#### 5. Secure copy files to your computer to view .html version fastqc files.
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
#### 1. Make a new skewer directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10
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

#### 5. Run fastqc on trimmed .fastq files to generate .fastqc files in html format.
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
##### 6b. View files in browser to check quality, then return to the sphinx terminal.

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

#### 3. Return to the 3_bwa directory and symbolically link trimmed fastq files from the 2_skewer directory.
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
* the -t flag specified the number of threads for `bwa` to use.
* the -bSh flag in `samtools` converts large SAM files into smaller, binary BAM files with headers.
* the -@ flag designates threads for the job in `samtools`.
* the -m flag designates memory for the job in `samtools`.
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
#### 7. Run the script, which will produce sorted BAM files, using bash with the following command:
```
bash bwa.sh
```
#### 8. To assess mapped reads, run the following for loop to create mapping stats files with `samtools`.
```
nano samtools_flagstat.sh

for file in *_sorted.bam ; do basename=$(echo $file | sed 's/-trimmed_sorted.bam//') ; samtools flagstat $basename-trimmed_sorted.bam > $basename-trimmed_sorted.stats ; done
```
