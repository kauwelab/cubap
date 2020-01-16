import csv
import argparse
from Bio import SeqIO

def parseArgs():
    '''
	Argument parsing is done.
	Required to have an input file.
	'''
    parser = argparse.ArgumentParser(description='Count Codon Usage.')
    parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
    args = parser.parse_args()
    return args

def invalid_codon(codon):
    '''
    Assesses validity of codon.
    Returns True if codon is invalid.
    '''
    for c in codon:
        if c != 'A' and c != 'T' and c != 'G' and c != 'C':
            return True
        elif c == '\n':
            return True

def get_next_semicolon(description):
    '''
    Returns the index of the next ';' in the gene description
    '''
    try:
        return description.index(';')
    except ValueError:
        return description.index('\\')

def get_gene_name(description):
    '''
    Extracts the gene name and NP number from the description
    '''
    geneName = description[description.find("gene=") + 5 :
        description.find("gene=") + 5 +
        get_next_semicolon(description[description.find("gene=") + 5 : ])]
    try:
        np = ' ' + description[description.find("Name=") + 5 :
            description.find("Name=") + 5 +
            get_next_semicolon(description[description.find("Name=") + 5 : ])]
    except:
        np = ''
    return geneName, np

def get_longest_isoform(genesDict):
    for gene,isoforms in genesDict.items():
        longestIsoform = ''
        longestLength = 0
        for isoform,length in genesDict[gene].items():
            if length > longestLength:
                longestLength = length
                longestIsoform = isoform
        print(gene, longestIsoform)

def sample_count(sample, sampleName):
    '''
    Reads codon usage data of every gene into a dict
    '''
    for seq_record in SeqIO.parse(sample, "fasta"):
        seq = seq_record.seq
        description = seq_record.description
        if description.find('exception') == -1 and description.find('partial=true') == -1 and len(seq) % 3 == 0:
            geneName, isoformNumber = get_gene_name(description)
            if geneName not in genesDict.keys():
                genesDict[geneName] = isoformDict.copy()
            genesDict[geneName][isoformNumber] = len(seq)

if __name__ =='__main__':
    '''
    Main.
    '''
    genesDict = {}
    isoformDict = {}
    args = parseArgs()
    sampleFile = str(args.input[0])
    sampleName = sampleFile[sampleFile.find('/') + 1:]

    sample_count(sampleFile, sampleName)
    get_longest_isoform(genesDict)