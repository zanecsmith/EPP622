# Test 2


#### Make a directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/
` cd /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/ ` \
` mkdir zsmith10 `

---
### Fastq QC [Directory: 1_fastqc]
----
#### 1. Make a new fastqc directory.
` mkdir 1_fastqc ` \
` cd 1_fastqc `

#### 2. Symbolically link Solenopsis invicta sequence files to 1_fastqc directory.
` ln -s ../../../raw_data/solenopsis_invicta_test2/*fastq . `

#### 3. Load fastqc in Spack.
` spack load fastqc `

#### 4. Process .fastq files into .fastqc files.
` for file in *.fastq; do fastqc $file; done `

#### 5. Secure copy files to your computer to view .html version fastqc files.
##### 5a. Open a new local terminal & store files sensibly.
` cd desktop ` \
` mkdir test2_1_fastqc ` \
` cd test2_ 1_fastqc ` \
` scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/1_fastqc/*html' . `
##### 5b. View files in browser to check quality, then return to sphinx terminal.

---
### Skewer [Directory: 2_skewer]
---
#### 1. Make a new skewer directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10
` cd ../ ` \
` mkdir 2_skewer ` \
` cd 2_skewer `

#### 2. Symbolically link .fastq files in the 2_skewer file.
` ln -s ../../../raw_data/solenopsis_invicta_test2/*fastq . `

#### 3. Search for adapter sequence counts in each file using the grep command. The -c flag will print counts of occurences.
##### Note: Illumina adapter sequence: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
` grep -c 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG' *.fastq `
###### Adapters found in each respective file for this test: 3459, 1106, 3105, and 964.

#### 4. Trim all .fastq files of synthetic adapter sequences using a for loop.
` for f in *fastq; do BASE=$( basename $f | sed 's/_1.fastq//g'); echo $BASE;  /sphinx_local/software/skewer/skewer -t 2 -l 95 -Q 30 -x AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG $f -o $BASE; done `
##### *the BASE=$( basename $f | sed 's/_1.fastq//g') function will define the basename of a file and remove the _1.fastq from all file names. Skewer will automatically add -trimmed.fastq to the base file names after analysis.
##### **The -t flag defines the number of threads for Skewer to use
##### **The -l flag defines minimum sequence lengths to keep
##### ** The -Q flag defines the minimum Phred quality score for seqeuences

#### 5. Run fastqc on trimmed .fastq files to generate .fastqc files in html format.
fastqc *trimmed.fastq

#### 6. Secure copy trimmed .fastqc.html files to your home computer's terminal.
` scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/2_skewer/*html' . `
##### 5b. View files in browser to check quality, then return to sphinx terminal.
