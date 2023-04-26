
import argparse
import pysam
import vcf as pyvcf
from collections import namedtuple

# For printing some results of the selected sample calling
SAMPLE_ID = 756
print("========================= python script 1 start =========================")
# Utility functions
def get_pos_idx(pos, aligned_pairs):
    idx = 0
    for pos_on_read, pos_on_ref in aligned_pairs:
        if pos_on_ref is None:
            idx += 1
            continue
        if pos-1 == pos_on_ref:   # Then pos_on_ref exists; pos-1 because 1-based index
            return idx
        idx += 1

    # Then pos_on_ref doesn't exists (None)
    # print("The position for this calling can't be found on the reference sequence (None)")
    return -1

def add_newFields_to(old_callData, RRNs: str, ARNs: str):
    new_CallData = namedtuple('CallData',old_callData._fields + ('RRN','ARN'))
    return new_CallData(*old_callData, RRN=RRNs, ARN=ARNs)

# Commandline Parser
parser = argparse.ArgumentParser(description='Argument that tells the script which bam file to process')
parser.add_argument('-bamFileName', type=str, help='file name for bam file')
parser.add_argument('-vcfFileName', type=str, help='file name for vcf file')
parser.add_argument('-outFileName', type=str, help='file name for output')
bam_Filename = parser.parse_args().bamFileName
vcf_Filename = parser.parse_args().vcfFileName
out_Filename = parser.parse_args().outFileName

print("bam_Filename:", bam_Filename)
print("vcf_Filename:", vcf_Filename)
print("out_Filename:", out_Filename)

''' Using the pyvcf lib to read as a vcf reader object'''
# vcf_Filename = "/data1/scaphase/hic-dovetail/var_" + file_ID + "_fixed.vcf"   # pyvcf cannot open .gz file
vcf_reader = pyvcf.Reader(open(vcf_Filename, 'r'))


''' Add new Format informations to header '''
print("================= Origional formats in header: =================")
for FORMAT in (vcf_reader.formats.values()):
    print('  ',FORMAT)
# define 2 new FORMAT fields called RRN and ARN
RRN_format = pyvcf.parser._Format('RRN', '.', 'String', 'Read Names of reads supporting the Ref Allele')
ARN_format = pyvcf.parser._Format('ARN', '.', 'String', 'Read Names of reads supporting the Alt Allele')
# The vcf.Reader object returned by pyvcf is not strictly read-only, so I can modify (but can't modify the file)
vcf_reader.formats['RRN'] = RRN_format
vcf_reader.formats['ARN'] = ARN_format
print("========= New formats added to reader object's header: =========")
for FORMAT in (vcf_reader.formats.values()):
    print('  ',FORMAT)
print("================================================================")

''' Using the pyvcf lib to read a vcf writer object'''
new_vcf_Filename = out_Filename
vcf_writer = pyvcf.Writer(open(new_vcf_Filename, 'w'), vcf_reader)  
''' Using the AlignmentFile method to read as a bam file object'''
# bam_Filename = "/data1/scaphase/hic-dovetail/hic_reads_" + file_ID + "_mapped.bam"
bam = pysam.AlignmentFile(bam_Filename,"rb")

print("Casting vcf file \n  ", 
        vcf_Filename,"\nto \n  ",
        new_vcf_Filename, "\nbased on bam file\n  ",
        bam_Filename)

print("================================================================")
longCalling = 0
n_records_left = 0
homo_or_non_diploids = 0
lowq = 0

''' Write the new vcf file '''
# Iterate over each record in the VCF file and add the new FORMAT and samples to each one
i = 0   # just want to print the first record
for record in vcf_reader:
    i += 1

    ''' = Filter all the low Qual callings / homo or non-diploids variants = '''
    if record.QUAL <= 30:
        lowq += 1       # QUAL <= 30
        continue

    GT = record.samples[0]['GT']
    diploids = ('0/1', '0/2', '1/0', '1/2', '2/0', '2/1')
    if  GT not in diploids:
        homo_or_non_diploids += 1       # Is not diploid vaiants or is homo_or_non_diploids
        continue

    ''' = Filter all the long callings = '''
    long_calling = False
    # For now we only consider single based variant callings 
    if len(record.REF) > 1:
        long_calling = True
    else:
        if len(record.ALT) > 1:
            long_calling = True
        else:
            for alt in record.ALT:
                if len(alt) > 1:
                    long_calling = True
                    break

    if long_calling:
        longCalling += 1
        # print(record.POS)
        continue
    else:
        n_records_left += 1

    ''' = Add the RRN and ARN FORMAT field to the record = '''
    record.add_format('RRN')
    record.add_format('ARN')
    '''Or'''
    # record.FORMAT += ':RRN:ARN'

    ''' = fetch and find the supported alleles = '''
    chrom = record.CHROM
    pos = record.POS
    RRNs = []   # Read Names supporting Reference allele
    ARNs = []   # Read Names supporting Alternative allele

    if i==SAMPLE_ID:        # show the aligned reads of a sample record    
        print("================ Aligned Reads of sample record ================")

    for read in bam.fetch(chrom, pos-1, pos):
        read_name = "_".join(read.query_name.split(":"))    # The read names include colon, we need to avoid confusion

        # Format for aligned_pairs:
        #   [(pos_on_read, pos_on_ref), ...]

        # Use get_aligned_pairs() to determine which allele the read supports
        aligned_pairs = read.get_aligned_pairs()
        pos_idx = get_pos_idx(pos, aligned_pairs)
        # print(aligned_pairs)

        if pos_idx == -1:   # Then pos_on_ref doesn't exists (None)
            # print("The position for this calling can't be found on the reference sequence (None)")
            continue
            
        aligned_pair = aligned_pairs[pos_idx]
        pos_on_read = aligned_pair[0]

        if pos_on_read is None:     # Then inserts, deletions, skipping may present
            # print("The position on the reading for this calling is none, inserts, deletions, skipping may present")
            continue
        
        if i==SAMPLE_ID:    # show the aligned reads of a sample record    
            print("read name:", read_name)
            print("query_sequence:", read.query_sequence)
            print("aligned pair:", aligned_pair)
            print("position on read:", pos_on_read)
            print("Allele on read:", read.query_sequence[pos_on_read])
            print("Reference Allele, i.e. record.REF[0]:", record.REF[0])           # REF may be like "AAT"
            print("Alternative Allele, i.e. record.ALT[0]:", record.ALT[0], '\n')   # ALT may be like "AATTCC". We only care SNPs in the first iteration

        if record.REF[0] == read.query_sequence[pos_on_read]:
            RRNs.append(read_name)
        elif record.ALT[0] == read.query_sequence[pos_on_read]:
            ARNs.append(read_name)
        # else:
        #     ARNs.append(read_name)
    ''' = done fetch and find the supported alleles = '''


    ''' = write the found alleles into sample field = '''
    RRN_str = ';'.join(RRNs)
    ARN_str = ';'.join(ARNs)

    # assign the new CallData object to the data attribute of the Call object
    call = record.samples[0]
    new_call_data = add_newFields_to(call.data, RRN_str, ARN_str)
    call.data = new_call_data
    vcf_writer.write_record(record)

    # Print the found RRNs and ARNs for the sample record
    if i == SAMPLE_ID:
        print("======================== Sample record ========================")
        print("Sample:\t i= ", i)
        print("CHROM:\t", record.CHROM)
        print("POS:\t", record.POS)
        print("ID:\t", record.ID)
        print("REF:\t", record.REF)
        print("ALT:\t", record.ALT)
        print("QUAL:\t", record.QUAL)
        print("FILTER:\t", record.FILTER)
        print("INFO:\t", record.INFO)
        print("FORMAT:\t", record.FORMAT)
        print("Sample:\t", record.samples[0])
        print("Sample Genotype:\t", record.samples[0]['GT'])
        print()
        print("Found RRNs:\n ", call.data[-2])
        print("Found ARNs:\n ", call.data[-1])

    # # Early termination for programming
    # if i>=1000:
    #     break

print("============================= Done =============================")
print("number of low quality callings (Qual <= 30):", lowq)
print("number of homozygous / non-diploid variants \n \
(GT other than 0/1, 1/0, 0/2, 2/0, 1/2, 2/1):", homo_or_non_diploids)
print("number of long callings (not single base):", longCalling)
print("number of records left:", n_records_left)
print("number of all records:", i)


''' Done writing new vcf file '''
# vcf_reader.close()  # not necessary, 'Reader' object doesn't have this attribute
vcf_writer.close()
bam.close()


print("========================== script end ==========================")

