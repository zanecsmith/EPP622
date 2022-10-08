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

#### 2. Link Solenopsis invicta sequence files to 1_fastqc directory.
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
##### 5b. View files in browser, then return to sphinx terminal

---
### Skewer [Directory: 2_skewer]
---
#### Make a new skewer directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10
` cd ../ ` \
` mkdir 2_skewer ` \
` cd 2_skewer `
