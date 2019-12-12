##########################


Make Allele Population Key File
makeAllelePopKey.py
Created by: Matt Hodgman
Email: matthodg@byu.edu


##########################


Purpose: to make a file that links to the main codon usage bias file and stores the population information for each allele, i.e. how many individuals from each subpopulation have a given allele.


##########################


ARGUMENT OPTIONS:

	-h	--help		show this help message and exit
	-i 	INPUT		Input Codon Usage Bias CSV File
	-o	OUTPUT		Output Directory/File Name
	-k	--keyfile	a csv key for identifying the sample's population, by default is "sampleKey.csv"

##########################


REQUIREMENTS:

makeAllelesPopKey.py uses Python version 3.7

Python libraries that must be installed include:
1. csv
2. argparse
3. pandas

If any of those libraries is not currently in your Python Path, use the following command:
pip install --user [library_name]
to install the library to your path.


##########################


Input Files:
This program requires one codon usage bias CSV file.

Output Files:
Output CSV file will be written to the same directory as makeAllelePopKey.py unless an alternate directory is included in the output CSV file name argument.


##########################


USAGE:

The -i argument is required, followed by the input CSV file. The -o argument is also required, followed by the desired name of the output CSV file.

By default, input CSV files are expected to be in the format outputted by the script SortCommonAlleles.py. 

Example usage:
python3 makeAllelesPopKey.py -i cub_data.csv -o cub_pop_data.csv


##########################
