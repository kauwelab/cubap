##########################


Codon Pairing Counter
codonPairingCounter.py
Created by: Justin Miller and Matt Hodgman
Emails: jmiller@byu.edu, matthodg@byu.edu


##########################


Purpose: count the number of identical and synonymous codon pairs


##########################


ARGUMENT OPTIONS:

	-h	--help		show this help message and exit
	-i 	INPUT		Input Fasta File
	-o	OUTPUT		Output Directory
	-f	FOOTPRINT	Ribosome Footprint window size
	-c	CO_TRNA		Flag to use co-tRNA codon pairing
	-b	BOTH		Flag to use both identical and co-tRNA codon pairing
	-p	POPULATION	Population of Current Input Fasta File
	-k	KEYFILE		Sample to Population Key CSV File
	-l	codon_table	Allows the user to specify an alternative codon table (BioPython). Default: Standard table


##########################


REQUIREMENTS:

codonPairingCounter.py uses Python version 3.7

Python libraries that must be installed include:
1. sys
2. csv
3. argparse
4. re
5. Bio.SeqIO,Bio.Seq,Bio.Alphabet

If any of those libraries is not currently in your Python Path, use the following command:
pip install --user [library_name]
to install the library to your path.


##########################


Input Files:
This program requires a single fasta file

Output Files:
Output file will be written to the same directory as codonPairingCounter.py or the output directory if provided. The file will have the same name as the input file but with the extension '.csv', '_aa.csv', or '_both.csv' depending on whether identical, co_tRNA, or both 
was passed as arguments respectively.


##########################


USAGE:

The -i argument is required, followed by the input faster file. The -p or -k arguments are also required followed by their respective values in order to correctly identify the population to which the sample belongs

By default, DNA sequences are expected. 

By default, the ribosome footprint size is 9 codons

Example usage:
python3 codonPairingCounter.py -i samples/HG00096 -p GBR -o samples_out/


##########################
