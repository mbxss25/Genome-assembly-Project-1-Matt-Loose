#Sahil Shelote Github repository (Rotation 1, group 6)

#Introduction
This project is about the identification of the genomic changes in the sequence of micro organism strains. However, we developed new assemblies with sample using minimaap, unicycler.
The quality check of assemblies were done with help of BUSCO and QUAST.
For the visualization and annotation we used flye, prokka and genovi.

#Things to note
group 6 made the repository for the codes we utilized for the sample set of 6.
The codes are wriiten in the respository with file name "final"
All code were runned in sbatch script.

#Installation and tools
conda - bash <conda-installer-name>-latest-MacOSX-x86_64.sh (for mac).
nanoplot - conda install -c bioconda nanoplot
minimap2 - conda install -c bioconda minimap.
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/minimap2/README.html)
Unicycler - conda install -c bioconda unicycler
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/unicycler/README.html)
QUAST - conda install -c bioconda quast
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/quast/README.html)
BUSCO - conda install -c bioconda busco
https://busco.ezlab.org/busco_userguide.html
Prokka - conda install -c bioconda prokka
[https://github.com/robotoD/GenoVi](https://github.com/tseemann/prokka)
flye - conda install flye
https://github.com/fenderglass/Flye/blob/flye/docs/INSTALL.md
Genovi - conda install -c bioconda genovi
https://github.com/robotoD/GenoVi

#DATA
#The data is provided in HPC at-

/workhere/students_2023/Matt_resources

-subfolders 
short reads (illumina data)
Long reads (Nanoporedata)

#Results expected (examples of output files).
#Nanoplot
Nanoplot (NanoPlot-report.html
NanoPlot_20240203_0110.log
NanoStats.txt
Non_weightedHistogramReadlength.html
Non_weightedHistogramReadlength.png
Non_weightedLogTransformed_HistogramReadlength.html
Non_weightedLogTransformed_HistogramReadlength.png
WeightedHistogramReadlength.html
WeightedHistogramReadlength.png
WeightedLogTransformed_HistogramReadlength.html
WeightedLogTransformed_HistogramReadlength.png
Yield_By_Length.html
Yield_By_Length.png

#Unicycler
001_spades_graph_k015.gfa
001_spades_graph_k029.gfa
001_spades_graph_k041.gfa
001_spades_graph_k049.gfa
001_spades_graph_k057.gfa
001_spades_graph_k063.gfa
001_spades_graph_k069.gfa
001_spades_graph_k073.gfa
002_depth_filter.gfa
003_overlaps_removed.gfa
004_long_read_assembly.gfa
005_bridges_applied.gfa
006_final_clean.gfa
assembly.fasta
assembly.gfa
miniasm_assembly/
unicycler.log

#Minimap
minimap_long_pf.fasta
minimap_long_pf.gfa
minimap_long_pf.paf.gz

#QUAST
aligned_stats/
basic_stats/
contigs_reports/
genome_stats/
icarus.html
icarus_viewers/
quast.log
report.html
report.pdf
report.tex
report.tsv
report.txt
transposed_report.tex
transposed_report.tsv
transposed_report.txt

#Busco
logs/
prodigal_output/
run_archaea_odb10/
short_summary.specific.archaea_odb10.bsc_hyb_pf.json
short_summary.specific.archaea_odb10.bsc_hyb_pf.txt

#Genovi
genovi_10_result-contig_1.png
genovi_10_result-contig_1.svg
genovi_10_result-contig_2.png
genovi_10_result-contig_2.svg
genovi_10_result-contig_3.png
genovi_10_result-contig_3.svg
genovi_10_result-contig_4.png
genovi_10_result-contig_4.svg
genovi_10_result.png
genovi_10_result.svg
genovi_10_result_COG_Classification.csv
genovi_10_result_COG_Classification.csv_percentage.png
genovi_10_result_COG_Classification_percentages.csv
genovi_10_result_COG_Histogram.png
genovi_10_result_Gral_Stats.csv

#Prokka
prokka_b10.err
prokka_b10.faa
prokka_b10.ffn
prokka_b10.fna
prokka_b10.fsa
prokka_b10.gbk
prokka_b10.gff
prokka_b10.log
prokka_b10.sqn
prokka_b10.tbl
prokka_b10.tsv
prokka_b10.txt

#Flye
00-assembly/
flye.log
params.json

#To see the results code following on terminal

python3 -m http.server 12121

on web input ip of HPC, with server. e.g http://10.102.161.12:12121

