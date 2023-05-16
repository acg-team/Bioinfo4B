# Genomics practice:

You will the following commands by one by through copy-pasting. Make sure there is no auto-correction happening with the words, dots, spaces, signs. Whenever you are asked a yes or no question or something like y/n comes up type y and press enter to continue.

- Download the Files folder to your laptop from here: https://github.com/acg-team/Bioinfo4B
- Initiate Anaconda by running the program
- Open a terminal: For mac: type terminal in your search and click on it. For windows: open a terminal in Anaconda prompt. 
- Create a directory called Genomics in your Documents and move everything inside the Files folder you downloaded in there. Use the following code
```
			cd Documents
			mkdir Genomics
			cp ../Downloads/Files/* Genomics
			cd Genomics 
			ls
```
- First, set the environment that will have all the packages we need. Thisa will take about 15 minutes.
```
			conda env create -f environment.yaml
			conda activate Env_Bioinfo4B
```
- Meanwhile, get your short reads using the sra-toolkit. For this, you'll actually install the program yourself. Go to the website: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit  Follow the instructions to download and install the toolkit. Open a tab in your terminal, make sure that the environment creation code is still running. 

Table for who works on which reads:
Daniela: SRR1583053

Daniel:  SRR1583055

Eric:    SRR1583057

Jonathan:SRR1583059

Kevin:   SRR1583061

Mark:    SRR1583063

```
			fastq-dump -I --split-files <type here the ID of the person’s genome assigned to you>
			cp <your SRR ID>.fastq reads.fastq
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
			cat consensus_short.fasta | sed -e 's/chrM/<name>/g’ > <name>.fasta
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
