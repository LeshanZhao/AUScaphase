# AUScaphase
====== Auburn Scaphase 202304 ========
====== Author: Leshan Zhao ===========
====== Advisor: Prof. Haynes Heaton ==

****** Description:
*   Scaphase (i.e. Scaffold Haplotype Phasing) is a bash script that
* integrates a serious of processes to perform haplotype phasing. It
* generates the counts of alleles based on the input reference
* sequence (a .fa / .fasta file) and hi-C reads (two .fq/.fastq/.gz files).

------ Usage:
- input:  1. ref.fa
-         2. read1.fq / read1.fq.gz (.gz preferred)
-         3. read2.fq / read2.fq.gz (.gz preferred)
-         4. scaffolds.agq
- output: allele_counts.out
- Sample usage 1:
-     bash scaphase.sh ref.fa read1.fq.gz read2.fq.gz scaffolds.agp output.out
- Sample usage 2:
-     bash scaphase.sh \\
/data1/scaphase/hic-dovetail/fAciRut3.PAT.20210629.primary_contigs.fa \
/data1/scaphase/hic-dovetail/28791_7_R1.fq.gz \
/data1/scaphase/hic-dovetail/28791_7_R2.fq.gz \
/data1/scaphase/hic-dovetail/scaffolds_FINAL.agp \
./scaphase_output.out

~~~~~~ Prerequisites:
The following softwares in the environment are required:
> minimap2-2.24 (r1122) (or similar versions)
> samtools 1.15.1 (or similar versions)
> freebayes v1.3.6 (or similar versions) \033[33m 
> python 3.6 (or earlier) (python 3.7 is incompatible with pyvcf) \033[0m

The following libraries in the python or conda environment are required:
> pysam 0.15.3 (or similar versions)
> pyvcf 0.6.8 (or similar versions)

