# Genomics:

## Preparation 
Here are your ZHAW credentials. Connect to the ZHAW network, using these. These are also your HPC credentials. In your first time, you will be asked to change the password for the HPC.


Xavier:	delfoxav	

Hiu Man:	usrn0016	TweJ$3Xsrc

Lova: usrn0017	sW4$Hkcrmt

Rebeca:	usrn0018	py3mLhqx!U

Hana: zemlihan	ks#@8MK94U

Besmira:	 sabn	

Philipp:	usrn0019	ks#@8MK94U

Emmanuelle:	usrn0020	D8N7YsuhkO

![image](https://github.com/acg-team/Bioinfo4B/assets/26571015/ac0e4a70-0ba8-40f1-aa5e-6efd0b25273b)

First, open the Visual Studio Code. In your first time, you will choose the ssh remote connection option and this will set up your connection via installing some necassary stuff. Once connected, open a new terminal using the menu. In the terminal, write the following code to move to the course folder in the terminal. This is where you will be working.

```
cd /cfs/earth/scratch/icls/shared/bioinfo4beginners
```
Check out the folders here with the command ls.

```
ls
```
Move to the Genomics folder. 

```
cd Genomics
```
Then create your own directory in there with your name. Replace myname with yours before writing ;). This is where you will be saving all your files.
```
mkdir tugce
```
Let's check if you did create your own folder:
```
ls
```

## Practice 

### Environment
Ok, let's begin! 
First, let's set the environment that will have all the packages we need. 
```
conda activate Env_Bioinfo4B
```


This diagram gives you an idea, what genomcs deal with :


![Steps-in-next-generation-sequencing-A-Extracted-DNA-is-randomly-broken-into-1000-bp](https://github.com/genombilim/2023/assets/37342417/6b4693c3-77b5-46e3-b74d-467425c933f8)

### download the reads
Then, get your short reads using the sra-toolkit. You'll send this as a job to the server. This way, we avoid programming on the main node, which would make the server too busy. First open the download.sh and look inside:
```
cat download.sh
```
Now, copy this to your own folder, don't forget to replce your name with mine.
```
cp download.sh tugce/tugce_download.sh
```
Let's see:
```
ls
```

Ok, now, we will download the reads. Use the table below in order to find the ID of the read you are assigned to:
Table for who works on which reads:

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
vim tugce_download.sh
```
Inside, you can use the up,down, left, right arrows to move in the file.
Press i in order to insert text in the file.
Go to where you see an SRR ID and change it to yours.
Press esc tab and then type :x in order to leave.

Run the code: 
```
sbatch tugce_download.sh
```
Check if it runs/is finished running etc.:
```
squeue
```

Once the job is done, you will have two fastq files in your folder. 

One for the forward reads that ends with _1 and one for the reverse reads, ending with _2. For simplicity, we will work with forward reads. Please change the name of the file as below:
```
mv SRR.... reads.fastq
```

![paired-end1](https://github.com/genombilim/2023/assets/37342417/3a672293-bb62-41b7-a361-0877512b8519)

Let's check out the file.
```
head reads.fastq
```

![Screen-Shot-2018-01-07-at-3 40 32-PM-1024x354](https://github.com/genombilim/2023/assets/37342417/1a2bed3d-f76d-442d-b74d-bf32657b3c3b)


Below is the link to the sra-toolkit. Explore it, when you have time:

https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit  

Also, here are the HPC Cluster's wiki pages on submitting and running a job:

https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Getting_started:Submit_a_job

https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Workload_management:Workload_management_(Slurm):Running_a_job_on_the_cluster


### Filter and trim
 
We’ll trim and filter the reads using the filter.sh script. Please copy this script to your own folder, as you change the name. Then we'll look in it again and send it as a job. 
```
cp ../filter.sh tugce_filter.sh
```
```
cat tugce_filter.sh
```
```
sbatch tugce_filter.sh
```

If you have time:

- Type command -h to explore what each paramater does and take note in your worksheet.
-  Check out the fastx-toolkit content with a pair: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html List two commands that look potentially useful, explore those and share with others. 

### Check out the Quality
![Screen-Shot-2018-01-07-at-1 36 09-PM-1024x713](https://github.com/genombilim/2023/assets/37342417/05a343ee-eed5-472c-86c0-08c1afa838ae)
- Compare the read sequence before and after trimming and filtering
```
fastqc reads.fastq
fastqc reads_trimmed.fastq
fastqc reads_filtered.fastq
```
Download those by right click and then double click on them to open them online. What difference do you notice? Why do you think?

This is the website of the program we used: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  
### Alignment to the reference genome
![image003](https://github.com/genombilim/2023/assets/37342417/e78cd5cd-4c55-4a0c-ac11-53b8b45e6a6b)
- Download the reference genome for the human mitochondria. For this, go to UCSC genome browser, choose Download, Human, Chromosomes and find the mitochondrial genome sequence. Copy-paste below where you see the word link:
```
python -m wget ‘link’
gunzip chrM.fa.gz
cp chrM.fa ref.fasta
  ```    
- Align the reads to the reference using the script align.sh. As usual, copy it to your folder and run it.

- You will notice that there is an indexing step. This figure is a visual depiction of the indexing:
- <img width="578" alt="index_kmer" src="https://github.com/genombilim/2023/assets/37342417/aa5fae6f-a6b0-4cc0-a2fc-12feddf0c7f9">

- Visualize the aligned sequences
 ```   
samtools tview -d C alignment_sorted.bam ref.fasta
 ```

Here you may find more information on aligment files: https://samtools.github.io/hts-specs/SAMv1.pdf

And this is the website of the tools we use for alignment: http://samtools.sourceforge.net/samtools.shtml

### Visualize the sequence depth

First type the code below:
```
samtools depth alignment_sorted.bam > depth.csv
```
Download the file. Then open the script Plot Depth.ipynb in Google Colab by choosing the repository name as acg-team.

### Identify the Variants
We will use the variants.sh file for this. As usual, make your own copy and run it. Lasty, enter the code below in irder to check out the contents of the vcf file:
```
bcftools stats calls.vcf.gz
```     
Here is more information on the VCF files: http://samtools.github.io/hts-specs/VCFv4.2.pdf
