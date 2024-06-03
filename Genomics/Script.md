# Genomics:

This diagram gives you an idea, what genomics deals with :

![Steps-in-next-generation-sequencing-A-Extracted-DNA-is-randomly-broken-into-1000-bp](https://github.com/genombilim/2023/assets/37342417/6b4693c3-77b5-46e3-b74d-467425c933f8)

The high throughput of next-generation sequencing (NGS) methods generates such large amounts of data that handling/filtering this data has become a major challenge during studies. Researchers also need to use more powerful computers or even take it one step further and use computing clusters. These new demands have given rise to bioinformatics. A bioinformatician is a Jack of all trades and bridges the gap between biology and computer science. Many Universities now offer specialized bioinformatics undergraduate and graduate programs but there are just as many self-taught bioinformaticians out there that originally come from a biology-only background.

It can’t be put any other way but, as in many other computer based sciences, performance comes before comfort. Most applications and programs can only be used via text commands in a terminal and they lack graphical user interfaces (GUIs). For users who are used to neatly designed user interfaces in Windows or Mac OS, using the terminal is associated with a steep learning curve. Many struggle with this at the beginning but nobody who wishes to work with next-generation sequencing data will get around learning at least some Linux/Unix terminal basics. Once you have worked with it for a while, however, you will get used to it and realize that there are good reasons for abandoning the comfort of GUIs:

	Everything has to be programmed, especially GUIs. New analysis methods for next-generation sequencing data are being developed and published on a weekly basis. This process would be slowed down significantly if the developer had to program a graphical user interface for every analysis software they write. As a consequence, your analyses will always lag behind current developments if you depend on commercial/GUI software.

	Graphical user interfaces occupy memory and CPU resources just to make things look good. Both memory and CPU time are major bottlenecks for analyses.

	Many of the files you are working on are in binary format or are too big to be displayed so they are not human readable anyways.

	Analyses are much easier to automate in a terminal. It is possible to “pipe” the output from one program directly into another program line by line as the analysis progresses. The second analysis step in the pipeline can therefore start even before the first one is finished. Once this is set up, you can start your entire pipeline and analyze your dataset from start to finish with one single keystroke.


In this practical, we aim to walk you through the entire process from raw data to the list of variants.

## Preparation 
Here are your ZHAW credentials. Connect to the ZHAW network, using these. These are also your HPC credentials.


Xavier:	delfoxav	

Hiu Man:	usrn0016	TweJ$3Xsrc

Lova: usrn0017	sW4$Hkcrmt

Rebeca:	usrn0018	py3mLhqx!U

Hana: zemlihan	ks#@8MK94U

Besmira:	 sabn	

Philipp:	usrn0019	ks#@8MK94U

Emmanuelle:	usrn0020	D8N7YsuhkO

![image](https://github.com/acg-team/Bioinfo4B/assets/26571015/ac0e4a70-0ba8-40f1-aa5e-6efd0b25273b)

**1) If this is your first time:** Visit https://idm.hpc.zhaw.ch. Login there with your initial password. You will be asked to set a new password. When this is set, logout.

**2) Login to the login node login-rhel8.hpc.zhaw.ch using SSH:** For this, open the Visual Studio Code. In your first time, you will choose the ssh remote connection option and this will set up your connection via installing some necassary stuff. Write the prompter your login name followed by the node name: <login.name>@login-rhel8.hpc.zhaw.ch It will ask you to enter your password in the prompter, do that. This connect you to the server.

Where you enter is your own working directory. You will run your codes here. Copy all the scripts here to run, save your outputs here. We also have a group directory. This is where all the scripts are stored.

**3)Set the environment:**  In order to run packages and programs we installed for you, you need to work in the environment, we prepared for you. For this, you need to perform multiple steps:

Once connected, you can open a new terminal using the menu. In the terminal, write the following code in order to initialize the shell for conda environment:
```
/cfs/earth/scratch/icls/shared/bioinfo4beginners/Genomics/miniconda3/bin/conda init
``` 
Logout to make changes active. 

Login to the login node login-rhel8.hpc.zhaw.ch again using SSH

Activate the environment by typing the code below:
```
conda activate /cfs/earth/scratch/icls/shared/bioinfo4beginners/Genomics/Env_Genomics
```
Copy the environment file to yourself:
```
cp /cfs/earth/scratch/icls/shared/bioinfo4beginners/Genomics/Env_Genomics.yml .
```


## Practice 

### Environment
Ok, let's begin! 
First, get all the scripts you need from the group folder and copy them to yours. 
```
cp /cfs/earth/scratch/icls/shared/bioinfo4beginners/Genomics/scripts/* .
```

### Download the reads
The study we are going to get our data from is “Maternal age effect and severe germ-line bottleneck in the inheritance of human mitochondrial DNA”. The study has uploaded the short reads it produced in the Sequence Read Archive (SRA). https://www.ncbi.nlm.nih.gov/sra?term=SRP047378

SRA is a bioinformatics database that provides a public repository for DNA sequencing data, especially the "short reads" generated by High-throughput sequencing, which are typically less than 1,000 base pairs in length. The SRA not only provides a place where researchers can archive their short read data, but also enables them to quickly access known data and their associated experimental descriptions (metadata) with pin-point accuracy.


Then, get your short reads using the sra-toolkit. For this, you need this package in your environment, as well. This is the only one that did not work with biopython. That's why we install it seperately. In order to add this in to you path, run:
```
export PATH=$PATH:$PWD/sratoolkit.3.1.1-ubuntu64/bin
```
Now check if it works:
```
which fastq-dump
```
You'll send this as a job to the server. This way, we avoid programming on the main node, which would make the server too busy. First open the download.sh and look inside:
```
cat download.sh
```
Use the table below in order to find the ID of the read you are assigned to:


Xavier: SRR1583064

Hiu Man:SRR1583063

Lavo: SRR1583062

Rebeca: SRR1583061

Hana: SRR1583060

Besmira: SRR1583059

Phillipp:SRR1583058

Manu: SRR1583057

Next, you will edit the file to download the sequence assigned to you. For this we will use the vim tool:
```
vim download.sh
```
Inside, you can use the up,down, left, right arrows to move in the file.

Press i in order to insert text in the file.

Go to where you see an SRR ID and change it to yours.

Press esc tab and then type :x in order to leave.

In case you don't want to save your changes or made some mistakes let's say: Press esc tab and then type :q! in order to leave without saving.

As you noticed, we will store the reports to another folder. So please first create this folder:
```
mkdir Report
```

And then run the code: 
```
sbatch download.sh
```
Check if it runs/is finished running etc.:
```
squeue
```

Once the job is done, you will have two fastq files in your folder. This is because we are dealing here with paired-end reads. One for the forward reads that ends with _1 and one for the reverse reads, ending with _2. For simplicity, we will work with forward reads. Please change the name of the file as below:
```
mv SRR...._1.fastq reads.fastq
```

![paired-end1](https://github.com/genombilim/2023/assets/37342417/3a672293-bb62-41b7-a361-0877512b8519)

Let's check out the file. This command will display the first 10 lines of the fastq file without opening the entire file. Remember, our files might be too big to open the entire file in a normal text editor as these editors usually load the entire file into memory at once which isn’t always possible. head (and its counterpart tail) is therefore a good command to briefly check if the files look like they’re supposed to. One read in a fastq file consists of 4 lines and looks something like this:
```
head reads.fastq
```

![Screen-Shot-2018-01-07-at-3 40 32-PM-1024x354](https://github.com/genombilim/2023/assets/37342417/1a2bed3d-f76d-442d-b74d-bf32657b3c3b)
the first line is the name of the sequence, the second line depicts the actual reads, the third line gives a sequence identifier and the length of that read, and the fourth line encodes the Phred scaled quality scores. 


Also check the Report folder. What do you see?

Below is the link to the sra-toolkit. Explore it, when you have time:

https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit  

And, here are the HPC Cluster's wiki pages on submitting and running a job:

https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Getting_started:Submit_a_job

https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Workload_management:Workload_management_(Slurm):Running_a_job_on_the_cluster


### Trim and filter
 
The large variation in quality scores and the identified biases in sequence content per base should be dealt with before we proceed with our analyses. We will do this in two steps. First, we will remove the low quality bases at the beginning and the end of the reads (trimming). In a second filtering step we will then remove all reads that have a high number of low-quality base calls in the middle section. 

We’ll trim and filter the reads using the filter.sh script. 
```
cat trim_filter.sh
```
```
sbatch trim_filter.sh
```

For the trimming: You will have a new fastq file in your working directory. The command above means that we are removing the first 20 bases (-f 20) and that the last position to be kept in the read is 240 (-l 240).

The last command means that we are filtering the dataset for reads that have more than 95% (-p 95) of bases above a sequencing quality of 30 (-q 30). All reads that do not fulfil this criterion will be discarded.

If you have time explore what each paramater does and take note in your worksheet:
-  Check out the fastx-toolkit content with a pair: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html List two commands that look potentially useful, explore those and share with others. 

### Check out the Quality
![Screen-Shot-2018-01-07-at-1 36 09-PM-1024x713](https://github.com/genombilim/2023/assets/37342417/05a343ee-eed5-472c-86c0-08c1afa838ae)

Each character there represents a quality value. Phred quality scores Q are defined as a property, which is logarithmically related to the base-calling error probabilities P. For example, if Phred assigns a quality score of 30 to a base, the chances that this base is called incorrectly are 10^((-30)/10)=0.001. This means that 1 out of 1000 calls with this quality score is actually wrong. Most studies put the threshold of base quality at 20 or more ideally at 30.

- Compare the read sequence before and after trimming and filtering using the quality.sh script. Check the contents and run it.
- The script will create html files for each reads file you have. Download those by right click and then double click on them in your local devide to open them online. In these files, some general statistics about the sequences are displayed.
- The basic statistics are some very basic facts about your sequence library such as the number of reads, GC content and sequence lengths. The number of reads in the library is usually the first thing to look at. Although it is no guarantee for good data, a high read count is usually a good sign.
- Click on “Per base sequencing quality” on the left. A boxplot is displayed. On the x-axis you can see the position within the read and on the y-axis the phred scaled quality score distribution across all reads at this position.
•	What do you notice regarding the shape of the curve?

The pattern we see here is quite common for illumina short read libraries: The read quality is slightly lower at the beginning and at the end of the read. We see the consequences of these lower sequencing quality scores when we look at the “Per base sequence content” (select the report via the list on the left). This report illustrates the frequency of each base across all reads on this position.

 •	When looking at the report, what do you notice concerning the base frequencies and how does it connect to the overall quality graph above?

- What differences do you notice between plots? Why do you think?
  
•	We do not want to use bad quality base calls for our alignments and variant calling later on. Why? What would be the risk if we DID use all reads despite low sequencing quality?

•	Why do you think it makes sense to do the trimming first?
•	Could/should you trim more based on your graphs? What do you think?
This is the website of the program we used: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  
### Alignment to the reference genome
Before aligning the reads to the reference mitochondrial genome, we first need to download the human mitochondrial reference sequence. 

![image003](https://github.com/genombilim/2023/assets/37342417/e78cd5cd-4c55-4a0c-ac11-53b8b45e6a6b)
- Download the reference genome for the human mitochondria. For this, go to UCSC genome browser, choose Download, Human, Chromosomes and find the mitochondrial genome sequence. Copy-paste below where you see the word link:
```
python -m wget ‘link’
```
```
gunzip chrM.fa.gz
cp chrM.fa ref.fasta
```   
- Align the reads to the reference using the script align.sh.

- Aligning reads to a genome can be viewed in the general context as string matching. The goal is to find a pattern (the short-read) in a large text (the reference genome), allowing for mismatches and indels. Naively, you can scan the text for the pattern but this is inefficient. There are techniques to pre-process (or index) the text to make queries fast and also that can even compress the size of the text. Bowtie2 uses this technique for ultrafast and memory-efficient aligment of sequencing reads to long reference sequences, which also supports paired-end alignment modes. We use bowtie2 commands to align the short reads. 

First, you index your reference genome, so that reads are quickly aligned. 
- You will notice that there is an indexing step. This figure is a visual depiction of the indexing:
- <img width="578" alt="index_kmer" src="https://github.com/genombilim/2023/assets/37342417/aa5fae6f-a6b0-4cc0-a2fc-12feddf0c7f9">
Note that for each different toolbox like SAMTools or bowtie2, you should index your reference genome by using their own commands. 

Next you align your short reads to the reference genome, using the reference file and your filtered file. 
It will take some minutes to create the alignment file, which is a SAM file (Sequence Alignment/Map), a human-readable TAB-delimited text format consisting of a header and an alignment section. Each alignment line has 11 mandatory fields for essential alignment information such as mapping position, mapping quality, read quality etc. After completion of the alignment, a brief report will appear in the terminal. Check it to see the total alignment rate of your reads. Finally, you convert your sam file to a bam (Binary Alignment/Map) file, which is simply a binary version of the sam file.

- Visualize the aligned sequences
 ```   
samtools tview -d C alignment_sorted.bam ref.fasta
 ```
Here you may find more information on aligment files: https://samtools.github.io/hts-specs/SAMv1.pdf

And this is the website of the tools we use for alignment: http://samtools.sourceforge.net/samtools.shtml

### Visualize the sequence depth
Of course we are now interested in how much of the mitochondrial genome we actually cover with our reads (and how many times we cover each position in the genome). Follow the instructions below to create a scatterplot in python.
First type the code below:
```
samtools depth alignment_sorted.bam > depth.csv
```
Download the file. Then open the script Plot Depth.ipynb in Google Colab by choosing the repository name as acg-team.

•	What is the lowest depth in your dataset? What does it mean?

### Identify the Variants
Last step is mapping the variants between the reads and the reference genome. To this end, we will use SAMTools which provides various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. 
The aligner cannot always assign a read to its point of origin with high confidence. For instance, a read that originated inside a repeat element might align equally well to many occurrences of the element throughout the genome, leaving the aligner with no basis for preferring one over the others.

Aligners characterize their degree of confidence by reporting mapping quality, just like read quality scores. A mapping quality of 10 or less indicates that there is at least a 1 in 10 chance that the read truly originated elsewhere. Ideally, you will filter out bases with read qualities lower than 30 but also with mapping qualities lower than 10. In our case, the alignments are already good, so we will skip this step. 

We will use the variants.sh script for this. This will create a variant call file and some plots of the variants under a directory called Plots. Check the substitutions.0.png plot. Compare with each other. Would it be possible to find similar or different patterns among your samples?
Here is more information on the VCF files: http://samtools.github.io/hts-specs/VCFv4.2.pdf
