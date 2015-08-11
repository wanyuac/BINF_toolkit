"""
This script parses an ENA metadata file in XML format and prints a subset of information.

Usage: python parse_ENA_sampleInfo_XML.py ERP000909.xml > samples.txt

Input: an XML file exported for a list of ERS accession numbers from ENA using the REST URLs API. For example, one can download an XML file
	   for sample ERS086023 using http://www.ebi.ac.uk/ena/data/view/ERS086023&display=xml.

Output: a tab-delimited text file containing information retrieved from the XML file.
		study_accession, sample_accession, secondary_sample_accession, experiment_accession, run_accession, Isolate_ID, Host, Place_of_isolation, Year_of_isolation

Author of this version: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Edition history: 6-7, 11 August 2015

Licence: GNU GPL 2.1
"""

import sys
import xml.etree.ElementTree as xmlTree

def get_domains(sample):
	study = BioSample = ERS = experiment = run = isolate = strain = host = place = year = "NA"  # default value of all fields
	for domain in sample:
		if domain.tag == "IDENTIFIERS":
			BioSample, ERS = sample[0][1].text, sample[0][0].text  # <tag>text</tag>
		if domain.tag == "SAMPLE_LINKS":
			study = sample[4][0][0][1].text  # visit nested elements with indices
			experiment = sample[4][1][0][1].text
			run = sample[4][2][0][1].text
		if domain.tag == "SAMPLE_ATTRIBUTES":  # This domain may be variable in terms of attributes
			for attribute in domain:
				if attribute[0].text == "collection_date":
					year = attribute[1].text
				elif attribute[0].text == "isolate":
					isolate = attribute[1].text
				elif attribute[0].text == "specific_host":
					host = attribute[1].text
				elif attribute[0].text == "country":
					place = attribute[1].text
				elif attribute[0].text == "strain":
					strain = attribute[1].text
	return [study, BioSample, ERS, experiment, run, isolate, strain, host, place, year]

def main():
	file = sys.argv[1]
	xml = xmlTree.parse(file).getroot()  # parse an XML into a tree of elements
	
	# print the header line
	print "\t".join(["study_accession", "sample_accession", "secondary_sample_accession", "experiment_accession", "run_accession", "Isolate_ID", "Strain", "Host", "Place_of_isolation", "Year_of_isolation"])
	for sample in xml:
		print "\t".join(get_domains(sample))
	return

if __name__ == '__main__':
	main()
