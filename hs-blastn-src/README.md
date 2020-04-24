Indroduction
----------------------------

The original source codes (version 0.0.5) are found in directory `v0.0.5`.

HS-BLASTN is a Nucleotide-Nucleotide aligner that searches nucleotide queries against a large nucleotide database.
It aims to speedup the NCBI-BLASTN (the default task, megablast).

In this version of `HS-BLASTN`, we employ many ideas from `MECAT`. We no longer pursue produceing the same results as `BLASTN`, but finding the alignments having highest scores.

Available
------------------------------

HS-BLASTN is written in C and C++, open source and released under GPLv3.

Download and Install
------------------------------

```shell
 $ git clone https://github.com/chenying2016/queries.git
 $ cd hs-blastn-src
 $ make
 ```
After compilation, the executable file `hs-blastn` is found is `queries/Linux-amd64/bin·.

How to use
------------------------------

 The tested queries can be downloaded from 
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


2. Perform database search
    ```shell
    $ hs-blastn -num_threads 20 -keep_db HSQueriesLarge.fa hg38.fa > large_queries_map.sam
    ```
    or
    ```shell
    $ hs-blastn -num_threads 20 -block_size 500 HSQueriesSmall.fa hg38.fa > small_queries_map.sam
    ```
    
Input Format
---------------------------

The usage of `hs-blastn` is
```shell
hs-blastn [OPTIONS] queries database > output.sam
```
Both `queries` and `database` can be one of the following formats:
* `FASTA` format file.
* `FASTQ` format file.
* File list format. For example:
```shell
$ cat query_file_list.txt
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190404-1.fail.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190404-1.pass.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190405-2.fail.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190405-2.pass.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190429-2.fail.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190429-2.pass.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190430-1.fail.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHONT20190430-1.pass.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHOON20190419-1.fail.fastq.gz
/data1/cy/ontsv/HG002_NA24385_Son/HacCall.WHOON20190419-1.pass.fastq.gz
$ hs-blastn -num_threads 20 query_file_list.txt hg38.fa > map.sam
```

Output Format
---------------------------

 HS-BLASTN currently supports 3 kinds of output formats (argument to the -outfmt option): 
 
 1. format 6 (tabular output without comments), 
 2. format 7 (tabular output with comments). 
 3. format 17 (SAM)
 
The default output format of `HS-BLASTN` is `17` (SAM).
 
HS-BLASTN Workflow
--------------------------

* Create a directory named `hbndb`. The name of this directory can be specified by the `-db_dir` option.
* Split the `queries` into volumes `Q1`, `Q2`, ..., `Qn`. Each volume is about `4G` residues. The number of DNA residues in one volume is specified by the `-max_query_vol_res` option. All the data are stored in the `hbndb` directory with `query` prefix.
* Split the `database` into volumes `S1`, `S2`, ..., `Sm`. Each volume is about `4G` residues. The number of DNA residues in one volume is specified by the `-max_subject_vol_res` options. All the data are stored in the `hbndb` directory with `subject` prefix.
* Align each `Qi` to each `Sj`.
* Delete the `hbndb` directory. If you want to keep this directory after the search, add the `-keep_db` flag to the `hs-blastn` command.
 
Citation
--------------------------

  Ying Chen, Weicai Ye, Yongdong Zhang, Yuesheng Xu. (2015) High speed BLASTN: an accelerated MegaBLAST search tool. Nucleic Acids Res., 43(16):7762-8.
 
Contact
--------------------------

 chenying2016@gmail.com
