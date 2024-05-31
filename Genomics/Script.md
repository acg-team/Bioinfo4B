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
cd tugce/tugce_download.sh
ls
```
Next, you will see how to edit the file. For this we will use the vim tool:
```
vim tugce_download.sh
```
Inside, you can use the up,down, left, right arrows to move in the file.
Press i in order to insert text in the file.
Press esc tab and then type x: in order to leave.

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

Run the code: 
```
sbatch tugce_download.sh
```
Check if it runs/is finished running etc.:
```
squeue
```

Below is the link to the sr-toolkit. Explore it, when you have time:

https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit  

Also, here are the HPC Cluster's wiki pages on submitting and running a job:

https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Getting_started:Submit_a_job

https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Workload_management:Workload_management_(Slurm):Running_a_job_on_the_cluster


### Filter and trim
 
We’ll trim and filter the reads using the filter.sh script.
```
			fastx_trimmer -f 20 -l 240 -i reads.fastq -o reads_trimmed.fastq
			fastq_quality_filter -q 30 -p 95 -i reads_trimmed.fastq -o reads_filtered.fastq
  ```    
- Type command -h to explore what each paramater does and take note in your worksheet.
-  Check out the fastx-toolkit content with a pair: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html List two commands that look potentially useful, explore those and share with others. 
 
- Compare the read sequence before and after trimming and filtering
```
			fastqc
  ```
  
  # Alignment to the reference genome
- Download the reference genome for the human mitochondria. For this, go to UCSC genome browser, choose Download, Human, Chromosomes and find the mitochondrial genome sequence. Copy-paste below where you see the word link:
```
			python -m wget ‘link’
			gunzip chrM.fa.gz
			cp chrM.fa ref.fasta
  ```    
- Align the reads to the reference:
```
			bowtie2-build ref.fasta ref
			bowtie2 -x ref -q reads_filtered.fastq -S alignment.sam
			samtools faidx ref.fasta
			samtools view -bt ref.fasta.fai alignment.sam > alignment.bam
			samtools sort alignment.bam -o alignment_sorted.bam
			samtools index alignment_sorted.bam
 ```    
 
- Visualize the aligned sequences
 ```   
			samtools tview -d C alignment_sorted.bam ref.fasta
 ```   
- Visualize the sequence depth. First type the code below, then open a console in Jupyter and follow the steps in Plot Depth.ipynb page in the Scripts folder in the GitHub repo.
```
			samtools depth alignment_sorted.bam > depth.csv
 ```     
  # Identify the Variants
- Identify the variants:
```
			bcftools mpileup -f ref.fasta alignment_sorted.bam | bcftools call -mv -Ov -o calls.vcf
			bcftools mpileup -f ref.fasta alignment_sorted.bam | bcftools call -mv -Oz -o calls.vcf.gz
  ```    
- Generate your first consensus sequence. Name it after yourself by replacing the <name> with your own.
```
			bcftools index calls.vcf.gz
			bcftools consensus calls.vcf.gz -f ref.fasta -o consensus.fasta
			fasta_formatter -i consensus.fasta -w 70 -o consensus_short.fasta
			cat consensus_short.fasta | sed -e 's/chrM/<name>/g' > <name>.fasta
  ```     
  # Draw a pylogenetic Tree
- Create a genomes files including yours and the population one:
```
			cat <your name>.fasta genomes/* > genomes.fasta
 ```     
- Create a phylogenetic tree
```
			clustalo -i genomes.fasta --outfmt=phylip -o genomes.phy
			dnaml
			figtree
```
To understand how dnaml works, go to its webiste https://evolution.genetics.washington.edu/phylip/doc/dnaml.html read through its assumptions. With a partner, discuss for each assumption how realistic they are and share what you discussed.
