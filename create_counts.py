
print("=================================== script start ===================================")

import argparse
import pysam
import vcf as pyvcf
# import gzip
SAMPLE = 1

# Commandline Parser
parser = argparse.ArgumentParser(description='Argument that tells the script which files to process')
parser.add_argument('-vcfFileName', type=str, help='file name for vcf file')
parser.add_argument('-agpFileName', type=str, help='file name for agp file')
vcf_Filename = parser.parse_args().vcfFileName
agp_Filename = parser.parse_args().agpFileName
print("vcf_Filename:", vcf_Filename)
print("agp_Filename:", agp_Filename)

''' Using the TabixFile method to read as a tbx file object'''
# Open the tabix indexed VCF file
tbi_Filename = vcf_Filename + ".tbi"
vcf_tbx = pysam.TabixFile(vcf_Filename, index=tbi_Filename)

''' Open the abp file and get the adjacent contigs within scaffolds '''
def get_scaffolds_from_agp(agp_Filename):
    agp = open(agp_Filename,'r')
    records = [line.strip('\n').split('\t') for line in agp.readlines()]
    scaffolds = {}
    for record in records:
        if record[0] not in scaffolds:
            scaffolds[record[0]] = []   # record[0] is scaffold name
        if record[4] == 'N':
            continue
        elif record[4] == 'W':
            contig_name = record[5]     # record[5] is contig name
        else:
            print("Something unexpected - Not N or W:", record[4])
            continue
        scaffolds[record[0]].append(contig_name)
    agp.close()   # should put this in the end
    return scaffolds

scaffolds = get_scaffolds_from_agp(agp_Filename)

# print(scaffolds)
print("======================================================================================")
print("================= Created hashmap 'scaffolds'                        =================")
print("================= Format:                                            =================")
print("=================    {'scaffold_1': ['contig1', 'contig2', ...],     =================")
print("=================     'scaffold_2': [...],                           =================")
print("=================     ...}                                           =================")
print("======================================================================================")
print("================= Sample (show only 5 scaffolds):                    =================")
sample_counter = 0
print("scaffolds:")
print("{")
for key, values in scaffolds.items():
    sample_counter += 1
    if sample_counter == 6:
        break
    print("'",key, "':")
    print(values,',')
print("...}")
print("======================================================================================")




print("\n================= Start creating hashmap all_contigs_readnames: =================\n")
sample_counter = 0
early_break = False


all_contigs_readnames = {}
for scaffold, contigs in scaffolds.items():
    if len(contigs) > 1:
        for contig_idx in range(len(contigs)):
            this_contig = contigs[contig_idx]

            ''' Fetch the call data for a specific contig '''
            readnames_in_this_contig = dict()
            # For each record, i.e., a new calling site
            for record in vcf_tbx.fetch(this_contig, parser=pysam.asTuple()):
                samples = record[-1]    # Tabix parser asVCF doesn't have field for samples, which is wierd
                RRNs = samples.split(':')[-2].split(';')
                ARNs = samples.split(':')[-1].split(';')
                for RRN in RRNs:
                    if RRN == '':
                        break
                    if RRN in readnames_in_this_contig:
                        readnames_in_this_contig[RRN].append('R')
                    else:
                        readnames_in_this_contig[RRN] = ['R']
                for ARN in ARNs:
                    if ARN == '':
                        break
                    if ARN in readnames_in_this_contig:
                        readnames_in_this_contig[ARN].append('A')
                    else:
                        readnames_in_this_contig[ARN] = ['A']

                ''' # --- Origional Implementation ---   
                for RRN in RRNs:
                    if RRN == '':
                        break
                    if RRN in readnames_in_this_contig:
                        if 'R' in readnames_in_this_contig[RRN]:
                            continue
                        else:
                            readnames_in_this_contig[RRN] = 'RA'
                    else:
                        readnames_in_this_contig[RRN] = 'R'
                for ARN in ARNs:
                    if ARN == '':
                        break
                    if ARN in readnames_in_this_contig:
                        if 'A' in readnames_in_this_contig[ARN]:
                            continue
                        else:
                            readnames_in_this_contig[ARN] = 'RA'
                    else:
                        readnames_in_this_contig[ARN] = 'A'
                    # --- Origional Implementation ---   
                '''
                # print("readnames_in_this_contig:")
                # print(readnames_in_this_contig)     # do something to the readnames_in_this_contig? like write it to a file or save it in a bigger dict

                sample_counter += 1
                if sample_counter == 1:
                    print('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv Hashmap (dictionary): all_contigs_readnames                       vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv Format:                                                           vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv    {"contig1": {"readname1": ["alleles1", "allele2", ...],        vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv                 "readname2": [...],                               vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv                 ...}                                              vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv     "contig2": {...}                                              vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv    ...}                                                           vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv')  
                    print('vvvvvvvvvvvvvvvvvvvvvv Sample: contig1: {readname: [alleles]}                            vvvvvvvvvvvvvvvvvvvvvv') 
                    print('vvvvvvvvvvvvvvvvvvvvvv                            found in first non-empty record of     vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvv                        the first contig of the first scaffold     vvvvvvvvvvvvvvvvvvvvvv')
                    print('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv')
                    print('all_contigs_readnames:')
                    print("{'",this_contig,"':")
                    print(readnames_in_this_contig)
                    print('}')
                    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^ Hashmap (dictionary): all_contigs_readnames                       ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^ Format:                                                           ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^    {"contig1": {"readname1": ["alleles1", "allele2", ...],        ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^                 "readname2": [...],                               ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^                 ...}                                              ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^     "contig2": {...}                                              ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^    ...}                                                           ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^ Sample: contig1: {readname: [alleles]}                            ^^^^^^^^^^^^^^^^^^^^^^') 
                    print('^^^^^^^^^^^^^^^^^^^^^^                            found in first non-empty record of     ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^                        the first contig of the first scaffold     ^^^^^^^^^^^^^^^^^^^^^^')
                    print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
                    print("\n================================ Still creating hashmap all_contigs_readnames: ================================\n")

                # if sample_counter == 100000:
                #     early_break = True
                if early_break:
                    break

            all_contigs_readnames[this_contig] = readnames_in_this_contig    #   save it in a bigger dict
            # hash_file.write(this_contig +'\n')
            # hash_file.write(str(readnames_in_this_contig))


            if early_break:
                break
        if early_break:
            break

print("\n========================== Hashmap done. Building counts: ==========================\n")


sample_counter = 0               
early_break = False
all_contigs_sites_count = dict()    # Counts_of_readnames_at_all_sites_of_all_contigs

for scaffold, contigs in scaffolds.items():

    if len(contigs) > 1:                                #   we only count for scaffolds with >= 2 contigs

        for contig_idx in range(len(contigs)-1):

            this_contig = contigs[contig_idx]           #   i.e., min_name
            next_contig = contigs[contig_idx+1]         #   i.e., max_name
            contig_pair = (this_contig, next_contig)
            all_contigs_sites_count[contig_pair] = [0, 0, 0, 0]     # initialization: [cis1, cis2, trans1, trans2]
            this_contig_pair_counts = all_contigs_sites_count[contig_pair]
            contig1_reads = all_contigs_readnames[this_contig]
            contig2_reads = all_contigs_readnames[next_contig]

            for readname, variants_ctg1 in contig1_reads.items():
                if readname in contig2_reads.keys():
                    variants_ctg2 = contig2_reads[readname]
                    
                    for allele1 in variants_ctg1:
                        for allele2 in variants_ctg2:
                            if allele1 == 'R' and allele2 == 'R':         # cis1
                                this_contig_pair_counts[0] += 1
                            elif allele1 == 'A' and allele2 == 'R':       # trans2

                                this_contig_pair_counts[3] += 1
                            elif allele1 == 'A' and allele2 == 'A':       # cis2

                                this_contig_pair_counts[1] += 1
                            elif allele1 == 'R' and allele2 == 'A':       # trans1

                                this_contig_pair_counts[2] += 1
                            else:                                         # Something wrong
                                print('allele1: ',allele1, 'allele2: ', allele2)
                                print('this_contig_pair_counts:', this_contig_pair_counts)
                                raise Exception("Something Unexpected")

            if early_break:
                break
        if early_break:
            break


print("\n================================ Done counting. ================================\n")
print("-------------------------------------------------------------------------------------------------------")
print("------------------ Hashmap (dictionary): all_contigs_sites_count                     ------------------")
print("------------------ Format:                                                           ------------------")
print("------------------       {(contig1, contig2): [n_cis1, n_cis2, n_trans1, n_trans2],  ------------------")
print("------------------        (contig2, contig3): [...],                                 ------------------")
print("------------------        ...}                                                       ------------------")
print("------------------ List every pair (min_contig_name, max_contig_name) below: v       ------------------")
print("-------------------------------------------------------------------------------------------------------")

print("all_contigs_sites_count = {")
for ctg_pair, counts in all_contigs_sites_count.items():
    print('\t contig pair: ', ctg_pair, "\t counts: ", counts)
print("}")

print("-------------------------------------------------------------------------------------------------------")
print("------------------ Listed every pair (min_contig_name, max_contig_name) above: ^     ------------------")
print("------------------ Hashmap (dictionary): all_contigs_sites_count                     ------------------")
print("------------------ Format:                                                           ------------------")
print("------------------       {(contig1, contig2): [n_cis1, n_cis2, n_trans1, n_trans2],  ------------------")
print("------------------        (contig2, contig3): [...],                                 ------------------")
print("------------------        ...}                                                       ------------------")
print("-------------------------------------------------------------------------------------------------------")
print('---- File ID (7 or 8, default is 7):', file_ID, ' ----\n')


vcf_tbx.close()     # should put this in the end

print("==================================== script end ====================================")




''' Create the new file to save the hash maps'''
# hash_Filename = "hashmap_" + file_ID + ".hmp"
# hash_file = open(hash_Filename,'w')    # or 'a' mode (append vs write)
# hash_file.close()   # should put this in the end
