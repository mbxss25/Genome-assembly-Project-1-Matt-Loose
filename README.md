# Project-1

Introduction: Group generated the strains of micro-organisms, however various genetic changes to the genomes were done but location of the changes were unknown. The long read data was mostly present. In order to know the changes assemblies were required. This respository consist of the codes for generating the assemblies from the long and short reads, along with the quality analysis on them.

Code are for generating:
-Short read assembly
-long read assembly
-A hybrid short and long read assembly

Following codes for Quality control of assemblies-
BUSCO
GENOVI
#————————————————————————————————————————--------------------------------------------------------------------------------#
This repository is for the samples allocated to group 6. i.e Group 6-sample 6 - illuminate S6 - Nanopore Barcode 06.

The data is provided in HPC at-

/workhere/students_2023/Matt_resources

-subfolders 
short reads (illumina data)
Long reads (Nanoporedata)
#——————————————————————————————————————————————————-------------------------------------------------------------------------#
#installation and activation of nanoplot in conda enviroment
conda create –n nanoplot_test –c bioconda nanoplot

#!/bin/bash


#SBATCH --job-name=grp6_nanoplot
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/mbxss25/slurm-%x-%j.out  #with your username

source $HOME/.bash_profile

#activate conda environment
conda activate /shared/conda/shared

# your command
NanoPlot


conda deactivate

#The expected result is creating the conda environment, installation and activation of the nanoplot within the conda environment.
#-------------------------------------------------------------------------------------------------------------------------------#

#Using nanoplot for checking the samples provided. The following code is for the illumina short reads provided in the dataset. The samples are differentiated between R1 and R2 as Forward and reverse. 

#nanoplot short reads (R1 and R2) 

#!/bin/bash

#SBATCH --job-name=short_reads_nanoplot
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out    # Replace 'payya4' with your username
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err     # Replace 'payya4' with your username

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run nanoplot on short reads (R1)
NanoPlot --fastq /workhere/students_2023/group6/illumina_short/R1_files/*.fastq.gz \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/nano_qc_R1_short

#run nanoplot on short reads (R2)
NanoPlot --fastq /workhere/students_2023/group6/illumina_short/R2_files/*.fastq.gz \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/nano_qc_R2_short

#deactivate conda env
conda deactivate

#The expected result running this code is 
#----------------------------------------------------------------------------------------#

#pass long reads

#!/bin/bash


#SBATCH --job-name=fastq_pass_long
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out	   #enter your username here
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err      #enter your username here

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run nanoplot on long reads (pass)
NanoPlot --fastq /workhere/students_2023/group6/nanopore_long/* \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/nano_qc_long_pass

#deactivte conda env
conda deactivate

#The expected result running this code is 
#----------------------------------------------------------------------------------------#

#fails long reads

#!/bin/bash

#SBATCH --job-name=fastq_fails_long
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out    #enter your username here
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err     #enter your username here

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run nanoplot on long reads
NanoPlot --fastq /workhere/students_2023/group6/fastq_fails_long/* \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/nano_qc_long_fails

#deactivate conda env
conda deactivate

#The expected result running this code is 
#----------------------------------------------------------------------------------------#

##Merging reads

#!/bin/bash

#SBATCH --job-name=merged_reads
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out    #enter your username here
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err



# Merge long reads (pass)
cat /workhere/students_2023/group6/nanopore_long/* > \
/workhere/students_2023/group6/merged_long_reads_pass.fastq.gz


# Merge long reads (fail)
cat /workhere/students_2023/group6/fastq_fails_long/* > \
/workhere/students_2023/group6/merged_long_reads_fail.fastq.gz


# Merge long reads (pass and fail)
cat /workhere/students_2023/group6/nanopore_long/* /workhere/students_2023/group6/fastq_fails_long/* > \
/workhere/students_2023/group6/merged_long_reads_pf.fastq.gz


# Merge short reads R1
cat /workhere/students_2023/group6/illumina_short/R1_files/* > \
/workhere/students_2023/group6/merged_R1_short.fastq.gz


# Merge short reads R2
cat /workhere/students_2023/group6/illumina_short/R2_files/* > \
/workhere/students_2023/group6/merged_R2_short.fastq.gz

#----------------------------------------------------------------------------------------#

##Minimap (for long assembly)

#long reads pass

#!/bin/bash

#SBATCH --job-name=minimap_long_assembly
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out     #enter your username here
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err


source $HOME/.bash_profile


# Activate conda environment
conda activate /shared/conda/shared

# Assembly with minimap2
minimap2 -t 8 -x ava-ont /workhere/students_2023/group6/merged_long_reads_pass.fastq.gz \
| gzip -1 > /workhere/students_2023/group6/minimap_long_pass/minimap_long_pass.paf

# Miniasm
#miniasm -f \
#/workhere/students_2023/group6/merged_long_reads_pass.fastq.gz \
#/workhere/students_2023/group6/minimap_long_pass/"long_assembly.paf.gz" INSERT > /workhere/students_2023/group6/minimap_long_pass/"long_assembly_pass.gfa"

# Convert the assembly to FASTA format
#awk '/^S/{print ">"$2"\n"$3}' /workhere/students_2023/group6/long_assembly.gfa > /workhere/students_2023/group6/long_assembly.fasta

# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#
#reads pass

#!/bin/bash

#SBATCH --job-name=minimap_long_pass
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err



source $HOME/.bash_profile


# Activate conda environment
conda activate /shared/conda/shared

# Assembly with minimap2
minimap2 -t 8 -x ava-ont /workhere/students_2023/group6/merged_long_reads_pass.fastq.gz /workhere/students_2023/group6/merged_long_reads_pass.fastq.gz | gzip -1 > /workhere/students_2023/group6/minimap_r$

# Miniasm
miniasm -f /workhere/students_2023/group6/merged_long_reads_pass.fastq.gz  /workhere/students_2023/group6/minimap_results/minimap_long_pass/minimap_long_pass.paf.gz > /workhere/students_2023/group6/minim$

# Convert the assembly to FASTA format
awk '/^S/{print ">"$2"\n"$3}' /workhere/students_2023/group6/minimap_results/minimap_long_pass/minimap_long_pass.gfa > /workhere/students_2023/group6/minimap_results/minimap_long_pass/minimap_long_pass.f$

# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#

#minimap pass and fail

#!/bin/bash

#SBATCH --job-name=minimap_long_passfail
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err




source $HOME/.bash_profile


# Activate conda environment
conda activate /shared/conda/shared

# Assembly with minimap2
minimap2 -t 8 -x ava-ont /workhere/students_2023/group6/merged_long_reads_pf.fastq.gz /workhere/students_2023/group6/merged_long_reads_pf.fastq.gz | gzip -1 > /workhere/students_2023/group6/minimap_resul$

# Miniasm
miniasm -f /workhere/students_2023/group6/merged_long_reads_pf.fastq.gz  /workhere/students_2023/group6/minimap_results/minimap_long_pf/minimap_long_pf.paf.gz > /workhere/students_2023/group6/minimap_res$

# Convert the assembly to FASTA format
awk '/^S/{print ">"$2"\n"$3}' /workhere/students_2023/group6/minimap_results/minimap_long_pf/minimap_long_pf.gfa > /workhere/students_2023/group6/minimap_results/minimap_long_pf/minimap_long_pf.fasta

# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#
#Unicycler short reads (merged R1 and R2)

#!/bin/bash

#SBATCH --job-name=unicycler_short
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err


# Activate conda environment
source $HOME/.bash_profile
conda activate /shared/conda/shared


#Assembly with unicycler (short)
unicycler -t 8 \
-1 /workhere/students_2023/group6/merged_R1_short.fastq.gz \
-2 /workhere/students_2023/group6/merged_R2_short.fastq.gz \
-o /workhere/students_2023/group6/unicycler_short_result


# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#

#unicycler long reads

#!/bin/bash

#SBATCH --job-name=unicycler_long
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err

# Activate conda environment
source $HOME/.bash_profile
conda activate /shared/conda/shared

# Assembly with Unicycler (long - pass)
unicycler -t 8 \
-l /workhere/students_2023/group6/merged_long_reads_pass.fastq.gz \
-o /workhere/students_2023/group6/unicycler_long_pass


# Deactivate conda environment
conda deactivate

#----------------------------------------------------------------------------------------#

##Unicycler (hybrid assembly)

#Unicycler hybrid pass

#!/bin/bash


#SBATCH --job-name=uni_hyb_p
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err


echo "the job name is: $JOB_NAME_ID"

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared


#run unicycler hybrid assembly (pass only)
unicycler -1 /workhere/students_2023/group6/merged_R1_short.fastq.gz \
-2 /workhere/students_2023/group6/merged_R2_short.fastq.gz \
-l /workhere/students_2023/group6/merged_long_reads_pass.fastq.gz \
-t 8 -o /workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pass


#check unicycler.log tail to see if assembly has completed
#deactivte conda env
conda deactivate
#----------------------------------------------------------------------------------------#

#Unicycler hybrid pass and fail

#!/bin/bash


#SBATCH --job-name=uni_hyb_pf_2
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=72:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err


echo "the job name is: $JOB_NAME_ID"

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run unicyler hybrid assembly (pass+fail)
unicycler -1 /workhere/students_2023/group6/merged_R1_short.fastq.gz \
-2 /workhere/students_2023/group6/merged_R2_short.fastq.gz \
-l /workhere/students_2023/group6/merged_long_reads_pf.fastq.gz \
-t 8 -o /workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pf

#check unicycler.log tail to see if assembly has completed
#deactivte conda env
conda deactivate

#----------------------------------------------------------------------------------------#
No code for bandage
#----------------------------------------------------------------------------------------#

##Quast (quality assesement of the assembly)

#!/bin/bash

#SBATCH --job-name=quast_analysis
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/payya4/slurm-quast-%x-%j.out
#SBATCH --error=/shared/home/payya4/slurm-quast-%x-%j.err


# Activate conda environment with QUAST installed
source $HOME/.bash_profile
conda activate /shared/conda/shared


#QUAST unicycler

# QUAST for short-read assembly
python /shared/conda/shared/bin/quast \
-o /workhere/students_2023/group6/quast_results/quast_short.unicycler \
-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
/workhere/students_2023/group6/unicycler_results/unicycler_short_result/assembly.fasta

# QUAST for long-read assembly (pass and fail)
python /shared/conda/shared/bin/quast \
-o /workhere/students_2023/group6/quast_results/quast_long_pf_unicycler \
-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
/workhere/students_2023/group6/unicycler_results/unicycler_long_result/assembly.fasta

# QUAST for long-read assembly (pass)
python /shared/conda/shared/bin/quast \
-o /workhere/students_2023/group6/quast_results/quast_long_pass_unicycler \
-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
/workhere/students_2023/group6/unicycler_results/unicycler_long_pass/assembly.fasta


# QUAST for hybrid assembly (short and long_pf)
#python /shared/conda/shared/bin/quast \
#-o /workhere/students_2023/group6/quast_results/quast_hybrid_pf \
#-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
#-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
#/workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pf/ADD MISSING PATHWAY

# QUAST for hybrid assembly (short and long_pass)
#python /shared/conda/shared/bin/quast \
#-o /workhere/students_2023/group6/quast_results/quast_hybrid_pass \
#-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
#-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
#/workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pass/ADD MISSING PATHWAY


#QUAST minimap

# QUAST for minimap long-read (pass and fail)
python /shared/conda/shared/bin/quast \
-o /workhere/students_2023/group6/quast_results/quast_minimap_pf \
-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
/workhere/students_2023/group6/minimap_results/minimap_long_pf/minimap_long_pf.fasta



# QUAST for minimap long-read (pass)
python /shared/conda/shared/bin/quast \
-o /workhere/students_2023/group6/quast_results/quast_minimap_pass \
-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
/workhere/students_2023/group6/minimap_results/minimap_long_pass/minimap_long_pass.fasta



# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#

##Quast final

#!/bin/bash


#SBATCH --job-name=quast_final_final
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-quast-%x-%j.out
#SBATCH --error=/shared/home/mbxjk6/slurm-quast-%x-%j.err


# Activate conda environment with QUAST installed
source $HOME/.bash_profile
conda activate /shared/conda/shared



# QUAST for all unicycler (short, long, hybrid) and minimap (long)
python /shared/conda/shared/bin/quast -o /workhere/students_2023/group6/quast_results/quast_final_final \
-r /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/GCF_033550185.1_ASM3355018v1_genomic.fna \
-g /workhere/students_2023/group6/ncbi_dataset/ncbi_dataset/data/GCF_033550185.1/genomic.gff \
/workhere/students_2023/group6/minimap_results/minimap_long_pass/minimap_long_pass.fasta \
/workhere/students_2023/group6/minimap_results/minimap_long_pf/minimap_long_pf.fasta \
/workhere/students_2023/group6/unicycler_results/unicycler_short_result/assembly.fasta \
/workhere/students_2023/group6/unicycler_results/unicycler_long_result/assembly.fasta \
/workhere/students_2023/group6/unicycler_results/unicycler_long_pass/assembly.fasta \
/workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pass/assembly.fasta \
/workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pf/assembly.fasta



# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#

##Busco (for analysis)

#hybrid assemblies 

#!/bin/bash

#SBATCH --job-name=busco_analysis_2
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/payya4/slurm-busco-%x-%j.out
#SBATCH --error=/shared/home/payya4/slurm-busco-%x-%j.err

# Activate conda environment with BUSCO installed
source $HOME/.bash_profile
conda activate /shared/conda/busco

# BUSCO for assembly hybrid pass and fail
busco \
-i /workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pf/assembly.fasta -o bsc_hyb_pf --out_path /workhere/students_2023/group6/busco_results_hyb -l archaea_odb10 -m genome

# BUSCO for assembly hybrid pass
busco \
-i /workhere/students_2023/group6/unicycler_results/unicycler_hybrid_pass/assembly.fasta -o bsc_hyb_pass --out_path /workhere/students_2023/group6/busco_results_hyb -l archaea_odb10 -m genome

# Deactivate conda environment
conda deactivate

#----------------------------------------------------------------------------------------#

##Busco analysis

#!/bin/bash

#SBATCH --job-name=busco_analysis
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/shared/home/payya4/slurm-busco-%x-%j.out
#SBATCH --error=/shared/home/payya4/slurm-busco-%x-%j.err

# Activate conda environment with BUSCO installed
source $HOME/.bash_profile
conda activate /shared/conda/busco

# BUSCO for short-read assembly
busco \
-i /workhere/students_2023/group6/unicycler_results/unicycler_short_result/assembly.fasta -o /workhere/students_2023/group6/busco_results/bsc_shortreads -l /shared/conda/busco_downloads/lineages/archaea_$


# BUSCO for long-read assembly (pass and fail)
busco \
-i /workhere/students_2023/group6/unicycler_results/unicycler_long_result/assembly.fasta -o /workhere/students_2023/group6/busco_results/bsc_longreads_pf -l /shared/conda/busco_downloads/lineages/archaea$

# BUSCO for long-read assembly (pass)
busco \
-i /workhere/students_2023/group6/unicycler_results/unicycler_long_pass/assembly.fasta -o /workhere/students_2023/group6/busco_results/bsc_longreads_pass -l /shared/conda/busco_downloads/lineages/archaea$

#BUSCO minimap

#BUSCO for minimap long-read (pass and fail)
busco \
-i /workhere/students_2023/group6/minimap_results/minimap_long_pf/minimap_long_pf.fasta -o /workhere/students_2023/group6/busco_results/mm_long_pf -l /shared/conda/busco_downloads/lineages/archaea_odb10 $

#BUSCO for minimap long-read (pass)
busco \
-i /workhere/students_2023/group6/minimap_results/minimap_long_pass/minimap_long_pass.fasta -o /workhere/students_2023/group6/busco_results/mm_long_pass -l /shared/conda/busco_downloads/lineages/archaea_$


# Deactivate conda environment
conda deactivate
#----------------------------------------------------------------------------------------#



#----------------------------------------------------------------------------------------#
#new reads

barcode 01 and 02
#!/bin/bash

#SBATCH --job-name=nanoplot_01_02
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err



source $HOME/.bash_profile


#activate conda env
conda activate /shared/conda/shared


#run nanoplot on barcode 01
NanoPlot --fastq /workhere/students_2023/group6/new_reads/raw_reads/barcode01/*fastq \
--threads 8 \
--plots kde dot \
-o /workhere/students_2023/group6/new_reads/nanoplot_results/barcode01_nanoplot


#run nanoplot on barcode02
NanoPlot --fastq /workhere/students_2023/group6/new_reads/raw_reads/barcode02_/* \
--threads 8 \
--plots kde dot \
-o /workhere/students_2023/group6/new_reads/nanoplot_results/barcode02_nanoplot


#deactivte conda
conda deactivate

#----------------------------------------------------------------------------------------#

barcode 06

#!/bin/bash


#SBATCH --job-name=nanoplot_b6
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out	  #enter your username here
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run nanoplot on long reads (pass)
NanoPlot --fastq /workhere/students_2023/group6/new_reads/raw_reads/barcode06/* \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/new_reads/nanoplot_results/barcode6_nanoplot

#deactivte conda env
conda deactivate

#----------------------------------------------------------------------------------------#
#barcode 09

#!/bin/bash


#SBATCH --job-name=nanoplot_b9
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out	  #enter your username here
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run nanoplot on long reads (pass)
NanoPlot --fastq /workhere/students_2023/group6/new_reads/raw_reads/barcode09/* \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/new_reads/nanoplot_results/barcode09_nanoplot

#deactivte conda env
conda deactivate

#----------------------------------------------------------------------------------------#
#barcode 10

#!/bin/bash


#SBATCH --job-name=nanoplot_b10
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/payya4/slurm-%x-%j.out	  #enter your username here
#SBATCH --error=/shared/home/payya4/slurm-%x-%j.err

source $HOME/.bash_profile

#activate conda env
conda activate /shared/conda/shared

#run nanoplot on long reads (pass)
NanoPlot --fastq /workhere/students_2023/group6/new_reads/raw_reads/barcode10/* \
         --threads 8 \
         --plots kde dot \
         -o /workhere/students_2023/group6/new_reads/nanoplot_results/barcode10_nanoplot

#deactivte conda env
conda deactivate

#----------------------------------------------------------------------------------------#
##Merging reads

#!/bin/bash

#SBATCH --job-name=merged_reads
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/mbxjk6/slurm-%x-%j.out    #enter your username here
#SBATCH --error=/shared/home/mbxjk6/slurm-%x-%j.err     #enter your username here


# Merge barcode01 reads 
cat /workhere/students_2023/group6/raw_reads/raw_reads/barcode01/* > \
/workhere/students_2023/group6/raw_reads/merged_barcode01.fastq.gz 

# Merge barcode02 reads 
cat /workhere/students_2023/group6/raw_reads/raw_reads/barcode02/* > \
/workhere/students_2023/group6/raw_reads/merged_barcode02.fastq.gz  


# Merge barcode06 reads 
cat /workhere/students_2023/group6/raw_reads/raw_reads/barcode06/* > \
/workhere/students_2023/group6/raw_reads/merged_barcode06_pass.fastq.gz


# Merge barcode09 reads 
cat /workhere/students_2023/group6/raw_reads/raw_reads/barcode09/* > \
/workhere/students_2023/group6/merged_barcode09_pass.fastq.gz


# Merge barcode10 
cat /workhere/students_2023/group6/raw_reads/raw_reads/barcode10/* /workhere/students_2023/group6/fastq_fails_long/* > \
/workhere/students_2023/group6/merged_barcode10_pass.fastq.gz


