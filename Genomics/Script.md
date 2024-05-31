# Genomics practice:

Here are your ZHAW credentials again. Connect to the ZHAW network, using these if you did not already. These are also your HPC credentials. In your first time, you will be asked to change the password for the HPC.


Xavier	Delfosse	delfoxav@students.zhaw.ch	delfoxav	

Hiu Man	Grisch	hiuman.grisch@yahoo.com 	usrn0016	TweJ$3Xsrc

Lovasoa Manuelle Sylviane	Rakotozafy	lovasoa.rakotozafy@uzh.ch	usrn0017	sW4$Hkcrmt

Rebeca	Scalco	rebeca.scalco@unibe.ch 	usrn0018	py3mLhqx!U

Hana	Zemlickova	zemlihan@students.zhaw.ch 	zemlihan	ks#@8MK94U

Besmira	Sabani	sabn@zhaw.ch	sabn	

Philipp	Heider	philipp-heider@outlook.com	usrn0019	ks#@8MK94U

Emmanuelle	Mohbat	manu.mohbat@gmail.com	usrn0020	D8N7YsuhkO

![image](https://github.com/acg-team/Bioinfo4B/assets/26571015/ac0e4a70-0ba8-40f1-aa5e-6efd0b25273b)



First, open the Visual Studio Code. In your first time, you will choose the ssh remote connection option and this will set up your connection via installing some necassary stuff. Once connected, a new terminal using the menu. In the terminal, write the following code to move to the course folder in the terminal.

```
cd /cfs/earth/scratch/icls/shared/bioinfo4beginners
```

- First, set the environment that will have all the packages we need. 
```
conda activate Env_Bioinfo4B
```
- Meanwhile, get your short reads using the sra-toolkit. For this, you'll actually install the program yourself. Go to the website: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit  Follow the instructions to download and install the toolkit. Open a tab in your terminal, make sure that the environment creation code is still running. 

Table for who works on which reads:

Xavier: SRR1583064
Hiu Man:SRR1583063
Lavo: SRR1583062
Rebeca: SRR1583061
Hana: SRR1583060
Besmira: SRR1583059
Phillipp:SRR1583058
Manu: SRR1583057

```
			fastq-dump -I --split-files <type here the ID of the person’s genome assigned to you>
			cp <your SRR ID>_1.fastq reads.fastq
  ```  
  # Filter the short reads from the sequencer
 
 - Check out the quality of your reads.
```
			fastqc
  ```
  
- We’ll trim and filter the reads. 
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
