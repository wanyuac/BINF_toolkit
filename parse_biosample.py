"""
Download records from the NCBI BioSample database and extract values of certain attributes. The main
output file is a tab-delimited file, which can be readily imported to Excel. The column order of this
file is determined by that of attribute names in the parameter '-a'.

Usage:
    samples='SAMN0001,SAMN0002,SAMN0003,SAMN0004,SAMN0005'
    attr='strain,collection_date,geo_loc_name,host,host_disease,isolation_source'
    python parse_biosample.py -i $samples -a $attr -e 'xx@xx.xx'
    python parse_biosample.py -i 'file:accessions.txt' -a $attr -e 'xx@xx.xx'  # One accession number a line in accession.txt

Copyright 2020 Yu Wan <wanyuac@126.com>
Publication: 15 March 2020
Licenced under GNU GPL 2.1.
"""

from __future__ import print_function
from Bio import Entrez
from argparse import ArgumentParser
import os
import sys
import time
import xml.etree.ElementTree as xmlTree


def parse_arguments():
    parser = ArgumentParser(description = 'Download BioSample records and extract specific attributes')
    parser.add_argument('-i', type = str, required = True, help = 'Comma-delimited accessions or a file of a list of acceessions')
    parser.add_argument('-a', type = str, required = True, help = 'Comma-delimited names of attributes whose values will be extracted')
    parser.add_argument('-e', type = str, required = True, help = 'User\'s email address required for accessing the NCBI database')
    parser.add_argument('-o', type = str, required = False, default = 'metadata.tsv', help = 'Filename for extracted attribute values')
    parser.add_argument('-d', type = str, required = False, default = './record', help = 'Directory name for downloaded BioSample records')
    parser.add_argument('-r', action = 'store_true', required = False, help = 'Flag it to override existing XML files')
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    check_outdir(args.d)
    attributes = args.a.split(',')
    
    # Download BioSample records
    print('Start to download BioSample records.')
    records = download_records(accessions = get_accession_numbers(args.i), email = args.e, \
                               out_dir = args.d, override = args.r)  # Create an XML file [acc].xml in the record directory
    
    # Parse the records
    print('\nStart to parse BioSample records.')
    f = open(args.o, 'w')  # Create the output file
    f.write('\t'.join(['BioSample'] + attributes) + '\n')  # Print the header
    for a, p in records.items():  # Accession number and XML path
        if p == None:
            continue
        else:
            extract_attributes(accession = a, xml_path = p, out_file = f, attrs = attributes)
    f.close()
    
    return


def extract_attributes(accession, xml_path, out_file, attrs):
    print("Parsing record %s:" % accession)
    xml = xmlTree.parse(xml_path).getroot()  # Read and parse the XML file
    parental_domain = xml[0]  # Tag: BioSample
    
    for d in parental_domain:
        if d.tag == 'Attributes':
            attr_dict = parse_attribute_domain(d)  # Write attributes into the output file
            parsed_attrs = attr_dict.keys()  # Note that not all attrs can always be found in parsed_attrs.
            new_line = [accession]
            for a in attrs:
                if a in parsed_attrs:
                    new_line.append(attr_dict[a])
                else:
                    print('    Warning: attribute %s is not found in record %s.' % (a, accession))
                    new_line.append('NA')  # A NULL space holder.
            out_file.write('\t'.join(new_line) + '\n')
            print('    Record %s has been successfully parsed.' % accession)
    
    return


def parse_attribute_domain(d):
    """
    Converts the 'Attributes' domain into a dictionary.
    """
    attrs = {}
    for a in d:
        # Compared to display_name and attribute_name, harmonized_name is expected to be conserved across records.
        attrs[a.attrib['harmonized_name']] = a.text  # For example: {isolation_source: urine}
    
    return attrs


def download_records(accessions, email, out_dir, override):
    """
    Download and save BioSample records as XML files, and return a dictionary {accession: file path}.
    Existing XML files are skipped by default, which saves lots of time. (Particularly when users re-run this script)
    """
    paths = {}
    Entrez.email = email
    
    for a in accessions:
        xml_path = os.path.join(out_dir, a + '.xml')
        if os.path.exists(xml_path):
            if override:
                print('Download and override existing record %s.' % a)
                try:
                    handle = Entrez.efetch(db = 'biosample', id = a, rettype = 'xml', retmode = 'text')  # In XML format
                    paths[a] = xml_path
                    f = open(xml_path, 'w')
                    f.write(handle.read())  # Save the current record as an XML file
                    f.close()
                except:
                    print('    Warning: record %s is not found. Skip this record.')
                    paths[a] = None
                time.sleep(1)
            else:
                print('Skip existing record %s.' % a)
                paths[a] = xml_path
        else:
            print('Download record %s.' % a)
            try:
                handle = Entrez.efetch(db = 'biosample', id = a, rettype = 'xml', retmode = 'text')  # In XML format
                paths[a] = xml_path
                f = open(xml_path, 'w')
                f.write(handle.read())  # Save the current record as an XML file
                f.close()
            except:
                print('    Warning: record %s is not found. Skip this record.')
                paths[a] = None
            time.sleep(1)
    
    return paths


def get_accession_numbers(s):
    """
    Returns a list of BioSample accession numbers
    """
    PREFIX = 'file:'
    if s.startswith(PREFIX):
        s = s[len(PREFIX) : ]  # Drop the prefix
        try:
            with open(s, 'r') as f:
                accs = f.read().splitlines()
        except:
            sys.exit('Error: accession file %s is not accessible.' % s)
    else:
        accs = s.split(',')
    print('    Altogether %i accession numbers have been imported.' % len(accs))

    return accs


def check_outdir(d):
    if os.path.exists(d):
        print('Output directory %s exists.' % d)
    else:
        os.system('mkdir ' + d)
        
    return


if __name__ == '__main__':
    main()
