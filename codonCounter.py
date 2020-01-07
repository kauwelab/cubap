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
    parser.add_argument('-o', '--output', type=str, help='the directory to which output files will be saved', required=False)
    parser.add_argument('-p', '--population', help='the population of the current sample being processed')
    parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population")
    parser.add_argument('--columnar', action='store_true', help='a flag to output results in a columnar csv', required=False)
    args = parser.parse_args()
    return args
    
def toCSV(sampleName, genesDict, genesDict2, outFile):
    '''
    writes the dict to a csv
    '''
    codons = ['AAA','AAC','AAG','AAT','ACA','ACC','ACG',
              'ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA',
              'CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT',
              'CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG',
              'GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC',
              'TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA',
              'TTC','TTG','TTT']
    header = ['Sample','Gene','Name'] + codons
    outFile = outFile + '.csv'
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for gene, codonDict in genesDict.items():
            codonRow = [sampleName, gene.split(' ')[0], gene.split(' ')[1]]
            for codon in codons:
                if codon in genesDict[gene].keys():
                    codonRow.append(genesDict[gene][codon])
                else:
                    codonRow.append(0)
            writer.writerow(codonRow)
            codonRow = []
        for gene, codonDict in genesDict2.items():
            codonRow = [sampleName, gene.split(' ')[0], gene.split(' ')[1]]
            for codon in codons:
                if codon in genesDict2[gene].keys():
                    codonRow.append(genesDict2[gene][codon])
                else:
                    codonRow.append(0)
            writer.writerow(codonRow)
            codonRow = []
    print(outFile,'created.')
    

def toColumnarCSV(sampleName, population, genesDict, genesDict2, outFile):
    '''
    writes the dict to a columanr csv
    '''
    codons = makeAllPossibleCodons()
    header = ['Sample', 'Population', 'Gene', 'Name', 'Codon', 'Count']
    outFile = outFile + '.csv'
    codonRow = []
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for gene, codonDict in genesDict.items():
            for codon in codons:
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                if codon in genesDict[gene].keys():
                    codonRow.append(codon)
                    codonRow.append(genesDict[gene][codon])
                else:
                    codonRow.append(codon)
                    codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
        for gene, codonDict in genesDict2.items():
            for codon in codons:
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                if codon in genesDict2[gene].keys():
                    codonRow.append(codon)
                    codonRow.append(genesDict2[gene][codon])
                else:
                    codonRow.append(codon)
                    codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
    print(outFile,'created.')

def makeDict():
    CodonsDict = {'Sample': '', 'Gene': ''}
    codons = makeAllPossibleCodons()
    for c in codons:
        CodonsDict[c] = 0
    return CodonsDict

def makeAllPossibleCodons():
	'''
	Returns a set of all 64 possible codons (DNA).
	'''
	from itertools import product
	codons = product("ACGT",repeat=3)
	codonsComb = set()
	for c in codons:
		codonsComb.add("".join(c))
	return codonsComb
            
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
    return geneName + np

def separateAlleles():
    '''
    Separates the alleles into separate dicts
    Results in two dicts, each containing an allele of every gene
    '''
    genesDictCopy = genesDict.copy()
    for gene in genesDictCopy:
        if "_A" in gene:
            genesDict2[gene[:-2]] = genesDict[gene].copy()
            genesDict2[gene[:-2]]['Gene'] = gene[:-2]
            del genesDict[gene]

def sample_count(sample, sampleName):
    '''
    Reads codon usage data of every gene into a dict
    '''
    CodonsDict['Sample'] = sampleName
    for seq_record in SeqIO.parse(sample, "fasta"):
        seq = seq_record.seq
        description = seq_record.description
        if description.find('exception') == -1 and description.find('partial=true') == -1 and len(seq) % 3 == 0:
            gene = get_gene_name(description)
            if gene in genesDict:
                gene = gene + "_A"
            CodonsDict['Gene'] = gene;
            genesDict[gene] = CodonsDict.copy()
            seq_size = len(seq) - (len(seq) % 3)
            index = 0
            while index < seq_size:
                codon = seq[index:index+3]
                while invalid_codon(codon):
                    index += 1
                    codon = seq[index:index+3]
                if codon in CodonsDict:
                    genesDict[gene][codon] += 1
                index += 3
    separateAlleles()

if __name__ =='__main__':
    '''
    Main.
    '''
    genesDict2 = {}
    genesDict = {}
    columns = ['Sample', 'Gene'] + list(makeAllPossibleCodons())
    CodonsDict = makeDict()
    args = parseArgs()
    sampleFile = str(args.input[0])
    sampleName = sampleFile[sampleFile.find('/') + 1:]
    if args.output is None:
        outFile = sampleName
    else:
        outFile = args.output + sampleName
    if args.population is None:
        import pandas as pd
        df_pop = pd.read_csv(args.keyFile)
        population = df_pop[df_pop['Sample'] == sampleName.rstrip()]['Subpopulation'].to_string().split()[1]
    else:
        population = args.population

    sample_count(sampleFile, sampleName)
    if args.columnar is True:
        toColumnarCSV(sampleName, population, genesDict, genesDict2, outFile)
    else:
        toCSV(sampleName, genesDict, genesDict2, outFile)
