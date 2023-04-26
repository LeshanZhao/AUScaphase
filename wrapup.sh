echo "
====== Auburn Scaphase 202304 ========
====== Author: Leshan Zhao ===========
====== Advisor: Prof. Haynes Heaton =="
echo "
****** Description:
*   Scaphase (i.e. Scaffold Haplotype Phasing) is a bash script that
* integrates a serious of processes to perform haplotype phasing. It
* generates the counts of alleles based on the input reference
* sequence (a .fa / .fasta file) and hi-C reads (two .fq/.fastq/.gz files)."
echo "
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
./scaphase_output.out"
           
echo -e "
~~~~~~ Prerequisites:
The following softwares in the environment are required:
> minimap2-2.24 (r1122) (or similar versions)
> samtools 1.15.1 (or similar versions)
> freebayes v1.3.6 (or similar versions) \033[33m 
> python 3.6 (or earlier) (python 3.7 is incompatible with pyvcf) \033[0m "
echo "
The following libraries in the python or conda environment are required:
> pysam 0.15.3 (or similar versions)
> pyvcf 0.6.8 (or similar versions)
" 
echo "
============== Start Script ==============
Starting: $(date)"

if [[ -z $1 ]];
then 
    echo "No parameter passed."
    echo "Four file: ref.fa read1.fq read2.fq scaffolds.agp needed."
    exit 1
else
    # echo "First Parameter = $1"

    if [[ -z $2 ]]
    then
        echo "No second parameter passed. read1.fq needed."
        exit 2
    else
        # echo "Second parameter = $2"

        if [[ -z $3 ]]
        then
            echo "No third parameter passed. read2.fq needed."
            exit 3
        else
            # echo "Third parameter = $3"

            if [[ -z $4 ]]
            then
                echo "No fourth parameter passed. scaffolds.agp needed."
                exit 4
            # else
                # echo "Fourth parameter = $4"
            fi
        fi
    fi
fi

ref_file=$1
echo ref_file: $ref_file

read_file_1=$2
echo read_file_1: $read_file_1

read_file_2=$3
echo read_file_2: $read_file_2

agp_file=$4
echo agp_file: $agp_file

echo "
~~~~~~ Step 1: minimap2 + samtools creating bamtemp.bam
Starting: $(date)"

minimap2 -ax sr -t 16 --no-pairing \
$ref_file \
$read_file_1 \
$read_file_2 \
| samtools view -bS - \
| samtools sort - \
> bamtemp.bam


echo "
~~~~~~ Step 2: freebayes creating vcftemp.vcf
Starting: $(date)"
freebayes -f \
$ref_file \
bamtemp.bam \
> vcftemp.vcf


echo "
~~~~~~ Step 3: python script creating var_with_readname_and_filtered_Qual+GT.vcf
Starting: $(date)"

python vcf_filter-QualGT+add-readnames.py \
-bamFileName vcftemp.vcf \
-vcfFileName bamtemp.bam \
-outFileName var_with_readname_and_filtered_Qual+GT.vcf \
> vcf_filter_log.log

echo "
~~~~~~ Step 4: python script creating counts allele_counts.out
Starting: $(date)"

bgzip -c var_with_readname_and_filtered_Qual+GT.vcf \
> var_with_readname_and_filtered_Qual+GT.vcf.gz

tabix -p vcf var_with_readname_and_filtered_Qual+GT.vcf.gz

python create_counts.py \
-vcfFileName var_with_readname_and_filtered_Qual+GT.vcf.gz \
-agpFileName $agp_file \
> allele_counts.out

# rm bamtemp.bam
# rm vcftemp.vcf


echo "
Finishing: $(date)
============== Script Done. =============="
exit 0
