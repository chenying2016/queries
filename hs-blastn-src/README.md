Indroduction
----------------------------

HS-BLASTN is a Nucleotide-Nucleotide aligner that searches nucleotide queries against a large nucleotide database.
It aims to speedup the NCBI-BLASTN (the default task, megablast) and produces the same results as NCBI-BLASTN.
HS-BLASTN first builds an FMD-index for the database (the index command), and then performs then MegaBLAST search algorithms (the align command).
HS-BLASTN adopts some source codes from BWA (FMD-index) and NCBI-BLAST (maskers, statistics, ungapped alignment algorithms and gapped alignment algorithms).
HS-BLASTN supports both FASTA and FASTQ format query files. The base qualities in FASTQ format files will not be considered.

Available
------------------------------

HS-BLASTN is written in C++, open source and released under GPLv3.

Download and Install
------------------------------

```shell
 git clone https://github.com/chenying2016/queries.git
 cd hs-blastn-src
 make
 ```
How to use
------------------------------

 In verion 0.0.5+ of HS-BLASTN, please use window masker version 2.3.0+ to generate the window masker file, otherwise the HS-BLASTN will not run.

 We take the experiemtns in the paper as an example to introduce how to use HS-BLASTN.
 The following 3 tools are provided by [NCBI-BLAST Version 2.3.0+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz): blastn, mekablastdb, windowmasker. 
 These three tools, together with hs-blastn, and all the files from directory queries, 
 should be put in the same directory. The tested queries can be downloaded from 
 https://github.com/chenying2016/experimental_queries. You should download it by klicking the "Download ZIP" button 
 at the right bottom corner. A file experimental_queries-master.zip is downloaded. Then you can decompress it using the 
 following shell commands:
 ```shell
 unzip experimental_queries-master.zip
 cd experimental_queries-master
 cat experimental_queries.tar.bz2.a* | tar xjv
 ```
 
1. Prepare the database
    Download the database human build 38 from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz.
    Unpack the file and put all the sequences in one file, hg38.fa, for example.
    
2. build a blastdb
 ```shell
    makeblastdb -in hg38.fa -input_type fasta -dbtype nucl
```

3. Generate a frequency counts file for the database using the WindowMasker
    ```shell
    windowmasker -in hg38.fa -infmt blastdb -mk_counts -out hg38.fa.counts
    ```

4. Convert the text format frequency file to obinary format
    ```shell
    windowmasker -in hg38.fa.counts -sformat obinary -out hg38.fa.counts.obinary -convert
    ```

5. Build an FMD-index
    ```shell
    hs-blastn index hg38.fa
    ```

6. Perform database search
    ```shell
    hs-blastn align -db hg38.fa -window_masker_db hg38.fa.counts.obinary -query query100.fa -out results_query_100.fa -outfmt 7
    ```
 
On Searching against the repeat-subregion-rich database
---------------------------

 When searching against a repeat-subregion-rich database, such as the human genomic database, 
 the WindowMasker is must be used for masking out these regions.
 Otherwise no performance advantages of HS-BLASTN is guaranteed.
 In fact, the existance of the repeated subsequences results in searching results that are consisted completely of these subsequences.
 Those results are not of biological interest.
 In addition, they will also make the searching time unexpected long.

Output results
---------------------------

 HS-BLASTN currently supports 3 kinds of output formats (argument to the -outfmt option): 
 
 1. format 0 (Standard pairwise alignments, the default output format), 
 2. format 6 (tabular output without comments), 
 3. format 7 (tabular output with comments). 
 
Citation
--------------------------

  Ying Chen, Weicai Ye, Yongdong Zhang, Yuesheng Xu. (2015) High speed BLASTN: an accelerated MegaBLAST search tool. Nucleic Acids Res., 43(16):7762-8.
 
Contact
--------------------------

 chenying2016@gmail.com
