from Bio import Entrez
from Bio import SeqIO
import textwrap
import time
import os
import argparse
import re
import csv
import sys
import Bio


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description = 'Pull taxa sequences: pull all sequences for genes of interest from specified species',\
            epilog = textwrap.dedent('''
            Additional information:
            This script takes a file of species (one per line) and genes of interest and retrieves them from NCBI.
            There are two main output files, a collection table and fasta files of the respective genes.
            The collection table lists all the retrieved sequences with the accession numbers. The
            collection names are formatted with the species name followed by their specimen voucher. This
            naming convention is retained throughout to ensure seamless concatenation.

            Several temporary files are also save to disk and deleted at the completion of the script. This
            can be disabled with the --verbose flag.
            
            pull_tax_species.py -i species_list.txt -e your.email@random.com -r "internal transcribed spacer",rpb1,rpb2'''))

parser.add_argument('-i', '--input', type = str, metavar = '', required = True, help = 'enter file of species names, one species per line')
parser.add_argument('-e', '--email', type = str, metavar = '', required = True, help = 'enter email address so NCBI can contact you if there is a problem')
parser.add_argument('-v', '--verbose', action = 'store_true', help = 'verbose, keep intermediate files')
parser.add_argument('-r', '--regions', nargs='+', required = True, help = 'regions to get, can be a comma seperated list. If regions contain spaces it must be in "" marks')

args = parser.parse_args()


Entrez.email = args.email

def parse_regions(regions):
    '''
    Parse the seperate regions from args.regions input. Returns a concatenated string
    with the regions in "" followed by [All regions] and OR/AND depending on the
    position in the string, e.g.
    "internal transcribed spacer"[All Fields] OR "matk"[All Fields] OR "rbcl"[All Fields] AND
    '''

    new_list = []
    cmd_str = ''
    # seperate the regions into new list
    for i in regions:
        newish_list=i.split(',')
        for j in newish_list:
            new_list.append(j)

    num_regions = len(new_list)
    count = 0
    out = ''
    # parse the regions into the correct format with [All fields] and AND/OR
    for region in new_list:
        count += 1
        if count < num_regions:
            out += '"{}"[All Fields] OR '.format(region)
        elif count == num_regions:
            out += '"{}"[All Fields] AND '.format(region)

    return out


def parse_species_file(infile):
    print('\nParsing species names')
    species = []
    infile_path = os.path.abspath(infile)
    with open(infile, 'r', encoding='utf-8-sig') as infile:
        if os.path.getsize(infile_path) == 0:
            print('File is empty!!!')
            sys.exit(1)
        duplicates = []
        for line in infile:
            line = line.strip().lower()
            if line not in species and len(line) > 0:
                species.append(line)
            elif line in species:
                print('# WARNING: {} is a duplicate'.format(line))
                if line not in duplicates:
                    duplicates.append(line)
        if duplicates == 0:
            print('No duplicates found\n')
        else:
            print('Parsing species names complete, {} duplicates found\n'.format(len(duplicates)))
    return species


def get_taxid(species):
    species_list = parse_species_file(species)
    txids = []
    missing_sp = 0
    print('Retreiving txids from NCBI')
    for species in species_list: ### Add try except for species that aren't present

        try:
            search = Entrez.esearch(term=species, db='taxonomy', retmode='xml')
            record = Entrez.read(search)

            if record['IdList'][0] not in txids:
                txids.append(record['IdList'][0])
        except IndexError:
            missing_sp += 1
            print('# WARNING: {} not at NCBI, is there a spelling mistake, synonym etc?'.format(species))

    print('Txids retreived. {} species were not available at NCBI.\n'.format(missing_sp))
    
    return txids


def pull_seqs_from_txid(sp):
    txids = get_taxid(sp)
    print('Retreiving sequences from NCBI')
    for txid in txids:
        regions = parse_regions(args.regions) + 'txid{}[Organism:exp] NOT "whole genome"[All fields] NOT mRNA[All fields] NOT "complete genome"[All fields] NOT "partial genome"[All Fields]'.format(txid)
        seq_search = Entrez.esearch(db='nuccore', term=regions, retmax=5000)
        seq_search_results = Entrez.read(seq_search)
        seq_search.close()
        uidList = ','.join(seq_search_results['IdList'])

        out_handle = open("Temp_File_taxa_hits_file.gb", "a")
        fetch_handle = Entrez.efetch(db="nucleotide", id=uidList, rettype="gbwithparts", retmode="text") #rettype="gb"
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    print('Sequences retreived from NCBI\n')


def parse_gb_full_files():
    '''
    This function reads in Temp_File_taxa_hits_file.gb and parses out the species, voucher
    sequence etc. It returns a dictionary only containing the species_col and the
    associated region:accessions plus the blast match info.
    '''
    print('Parsing GenBank files')
    infile = open("Temp_File_taxa_hits_file.gb")
    record_iterator = SeqIO.parse(infile, "gb")
    count = 0
    no_voucher_num = 0
    dic_gb = {}

    for record in record_iterator:
        try:
            count += 1
            organism = ''
            accession = record.name
            isolate = ''
            specimen_voucher = ''
            sequence = record.seq
            gene = ''
            header_gene = ''
            gene_for_table = ''
            header = record.description.lower()
            strain = ''
            clone = ''
            species_col = ''

            #Check if the sequence is ITS, LSU or SSU
            if 'internal transcribed spacer' in header or 'its1' in header or 'its2' in header or 'its region' in header:
                header_gene = 'ITS'
            elif 'external transcribed spacer' in header:
                header_gene = 'ETS'
            elif '28s' in header or 'large subunit' in header or '25s' in header and 'internal transcribed spacer' not in header:
                header_gene = 'LSU'
            elif '18s' in header or 'small subunit' in header and 'internal transcribed spacer' not in header and '5.8s' not in header:
                header_gene = 'SSU'
            elif 'psba-trnh' in header:
                header_gene = 'psbA-trnH'
            elif 'psbk-psbi' in header:
                header_gene = 'psbK-psbI'
            elif 'atpf-atph' in header:
                header_gene = 'atpF-atpH'
            elif 'trnl-trnf' in header:
                header_gene = 'trnL-trnF'
            elif 'atpb-rbcl' in header:
                header_gene = 'atpB-rbcL'
            else:
                header_gene = 'unknown'

            #Check if the sequence is from a type specimen, often contained in the header
            if 'type' in header:
                header_gene = header_gene + '-TYPE'

            #Often the ITS and LSU regions are concatenated together. Make a length cut off of
            #700 bp and try and account for some length variation by making the min len difference of 100 (800 total length)
            elif header_gene == 'ITS' and len(sequence) > 700 and len(sequence) - 700 >= 100:
                header_gene = 'ITS and LSU?'
            else:
                header_gene = header_gene

            #Try and catch cases where the collection is undescribed and also has the voucher in the
            #organism line e.g. Cortinarius sp. OTA00001
            if record.annotations['organism'].split(' ')[1] == 'sp.':
                organism = record.annotations['organism'].split(' ')
                organism = organism[:2]
                organism = '_'.join(organism).replace('.', '')
            #replace white space with '_'
            elif record.annotations['organism'].split(' ')[1] != 'sp.':
                organism = record.annotations['organism'].replace(' ', '_')


            #Find some sort of voucher and gene
            for feature in record.features:
                for key, value in feature.qualifiers.items():
                    if 'isolate' in key:
                        isolate = ''.join(value)
                    elif 'specimen_voucher' in key:
                        #[David] specimen voucher are coded as "XXXXX/XX" by convention, but
                        #some people put it as "XXXXX_XX". So, change the last "_" to "\".
                        specimen_voucher_temp = ''.join(value)
                        specimen_voucher = re.sub(r"(\S+)_(\S+)$", r"\1/\2", specimen_voucher_temp)
                    elif 'strain' in key:
                        strain = ''.join(value)
                    elif 'clone' in key:
                        clone = ''.join(value)
                    #Check what the gene is. If the sequence is from a single gene with multiple exons
                    #then only the gene name will taken. If there are multiple genes present in the sequence
                    #then the different gene names will be concatenated together e.g. trnL-trnF
                    elif 'gene' in key and len(gene) == 0:
                        gene = ''.join(value)
                    elif 'gene' in key and ''.join(value) == gene:
                        gene = gene
                    elif 'gene' in key and str(gene).find(''.join(value)) == -1:#len(gene) != 0:
                        gene = gene + str('-' + ''.join(value))

            #Check what vouchers are present and concatenate organism name and voucher together
            #the specimen voucher has precedence over isolate, strain, clone
            if len(specimen_voucher) != 0:
                species_col = organism + '__' + specimen_voucher
            elif len(specimen_voucher) == 0 and len(isolate) != 0:
                species_col = organism + '__' + isolate
            elif len(specimen_voucher) == 0 and len(isolate) == 0 and len(strain) != 0:
                species_col = organism + '__' + strain
            elif len(specimen_voucher) == 0 and len(isolate) == 0 and len(strain) == 0 and len(clone) != 0:
                species_col = organism + '__' + clone
            elif len(specimen_voucher) == 0 and len(isolate) == 0 and len(strain) == 0 and len(clone) == 0:
                no_voucher_num += 1
                species_col = organism + '__no_voucher_{}'.format(str(no_voucher_num))

            #replace all spaces in collection names, otherwise it will cause issues when the fasta file is imported into other software
            species_col = species_col.replace(' ', '-')

            #check that ITS not in header
            if len(gene) != 0: #and header_gene == 'unknown'
                gene_for_table = gene
            else:
                gene_for_table = header_gene

            #add species_col to dic_gb with the associated regions and accessions
            if species_col not in dic_gb:
                dic_gb[species_col] = {gene_for_table:[accession, str(sequence)]}
            elif species_col in dic_gb:
                dic_gb[species_col][gene_for_table] = [accession, str(sequence)]

        except Bio.Seq.UndefinedSequenceError:
            print("{} doesn't include a sequence, it might be a master record".format(record.name))
            continue
    print('GenBank files parsed\n')
    return dic_gb



def write_table_csv_fasta():
    '''
    Take the dictionary from parse_gb_full_files and writes a csv file with the blast matches
    Also writes a fasta file with all the sequences from the blast results.
    '''
    print('Making collections table')
    #get all the regions from parse_gb_full_files dictionary for the csv header
    dict_data = parse_gb_full_files()
    temp_dict={}

    # First parse dict_data to extract only the info needed for the csv file
    for key,value in dict_data.items():
        for region in value:
            gene = region
            accession = value[region][0]
            #print(gene, accession)

            if key not in temp_dict:
                temp_dict[key] = {region:accession}
            else:
                temp_dict[key][region]=accession


    csv_fields = ['Species and collection']

    for species_col, regions in temp_dict.items():
        for region in regions.keys():
            if region not in csv_fields:
                csv_fields.append(region)

    with open('collections_table.csv', 'w') as file:
        writer = csv.DictWriter(file, fieldnames = csv_fields, extrasaction = 'ignore')
        writer.writeheader()
        for key,value in sorted(temp_dict.items()):
            row = {'Species and collection':key}
            row.update(value)
            writer.writerow(row)
    write_fasta()
    print('Collections table written to collections_table.csv')


def write_fasta():
    record_dict = parse_gb_full_files()
    for record, regions in record_dict.items():
        for region in regions:
            accession = regions[region][0]
            sequence = regions[region][1]

            with open('{}_sequences.fasta'.format(region),'a') as fastafile:
                fastafile.write('>{} {} {}\n{}\n'.format(record, region, accession, sequence))



def remove_temp_files():
    dir_name = os.getcwd()
    test = os.listdir(dir_name)
    for item in test:
        if item.startswith("Temp_File_"):
            os.remove(os.path.join(dir_name, item))


if __name__ == '__main__':
    if args.verbose:
        start = time.time()
        pull_seqs_from_txid(args.input)
        parse_gb_full_files()
        write_table_csv_fasta()
        finish = time.time()
        sec_diff = round(finish - start, 1)
        min_dif = round(sec_diff/60, 1)
        print('Time in seconds {}, or {} minutes\n'.format(sec_diff, min_dif))
    else:
        start = time.time()
        pull_seqs_from_txid(args.input)
        parse_gb_full_files()
        write_table_csv_fasta()
        remove_temp_files()
        finish = time.time()
        sec_diff = round(finish - start, 1)
        min_dif = round(sec_diff/60, 1)
        print('Time in seconds {}, or {} minutes\n'.format(sec_diff, min_dif))
