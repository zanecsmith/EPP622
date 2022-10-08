## Test 2


#### Make a directory in /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/
``` cd /pickett_shared/teaching/EPP622_Fall2022/analysis_test2/ ```
``` mkdir zsmith10 ```

### Fastq QC
#### Make a new fastq directory
``` mkdir 1_fastq ```
``` cd 1_fastqc ```

#### Link Solenopsis invicta sequence files
``` ln -s ../../../raw_data/solenopsis_invicta_test2/*fastq . ```

#### Load fastqc in Spack
spack load fastqc

#### Process fastq files into fastqc files
``` for file in *.fastq; do fastqc $file; done ```

#### Secure copy files to your computer to view .html version fastqc files
##### Open a new local terminal & store files sensibly.
``` cd desktop ```
``` mkdir test2_1_fastqc ```
``` cd test2_ 1_fastqc ```
``` scp 'zsmith10@sphinx.ag.utk.edu:/pickett_shared/teaching/EPP622_Fall2022/analysis_test2/zsmith10/1_fastqc/*html' . ```
