import csv
import argparse
from Bio import SeqIO

def parseArgs():
    '''
	Argument parsing is done.
	Required to have an input file.
	'''
    parser = argparse.ArgumentParser(description='calculate nucleotide composition.')
    parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
    parser.add_argument('-o', '--output', type=str, help='the directory to which output files will be saved', required=False)
    parser.add_argument('-p', '--population', help='the population of the current sample being processed')
    parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population")
    parser.add_argument('-c', '--columnar', action='store_true', help='a flag to output results in a columnar csv', required=False)
    args = parser.parse_args()
    return args
    
def toCSV(sampleName, genesDict, genesDict2, outFile):
    '''
    writes the dict to a csv
    '''
    bases = ['A', 'T', 'G', 'C']
    header = ['Sample','Gene','Name'] + bases
    outFile = outFile + '.csv'
    baseDict = []
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for gene, baseDict in genesDict.items():
            baseRow = [sampleName, gene.split(' ')[0], gene.split(' ')[1]]
            for base in bases:
                if base in genesDict[gene].keys():
                    baseRow.append(genesDict[gene][base])
                else:
                    baseRow.append(0)
            writer.writerow(baseRow)
            baseRow = []
        for gene, baseDict in genesDict2.items():
            baseRow = [sampleName, gene.split(' ')[0], gene.split(' ')[1]]
            for base in bases:
                if base in genesDict2[gene].keys():
                    baseRow.append(genesDict2[gene][base])
                else:
                    baseRow.append(0)
            writer.writerow(baseRow)
            baseRow = []
    print(outFile,'created.')
    

def toColumnarCSV(sampleName, population, genesDict, genesDict2, outFile):
    '''
    writes the dict to a columanr csv
    '''
    bases = ['A', 'T', 'G', 'C']
    header = ['Sample', 'Population', 'Gene', 'Name', 'Base', 'Count']
    outFile = outFile + '.csv'
    baseRow = []
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for gene, baseDict in genesDict.items():
            for base in bases:
                baseRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                if base in genesDict[gene].keys():
                    baseRow.append(base)
                    baseRow.append(genesDict[gene][base])
                else:
                    baseRow.append(base)
                    baseRow.append(0)
                writer.writerow(baseRow)
                baseRow = []
        for gene, baseDict in genesDict2.items():
            for base in bases:
                baseRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                if base in genesDict2[gene].keys():
                    baseRow.append(base)
                    baseRow.append(genesDict2[gene][base])
                else:
                    baseRow.append(base)
                    baseRow.append(0)
                writer.writerow(baseRow)
                baseRow = []
    print(outFile,'created.')

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
    BasesDict['Sample'] = sampleName
    for seq_record in SeqIO.parse(sample, "fasta"):
        seq = seq_record.seq
        description = seq_record.description
        if description.find('exception') == -1 and description.find('partial=true') == -1 and len(seq) % 3 == 0:
            gene = get_gene_name(description)
            if gene in genesDict:
                gene = gene + "_A"
            BasesDict['Gene'] = gene;
            genesDict[gene] = BasesDict.copy()
            for base in seq_record.seq:
                if base in BasesDict:
                    genesDict[gene][base] += 1
                else:
                    print("Invalid Base:", base, "Gene:", gene)
    separateAlleles()

if __name__ =='__main__':
    '''
    Main.
    '''
    genesDict2 = {}
    genesDict = {}
    columns = ['Sample', 'Gene', 'A', 'T', 'G', 'C']
    BasesDict = {'Sample': '', 'Gene': '', 'A': 0, 'T': 0, 'G': 0, 'C': 0}
    args = parseArgs()
    sampleFile = str(args.input[0])
    if sampleFile.find('.') == -1:
	sampleName = sampleFile[-7:]
    else:
	sampleName = sampleFile[sampleFile.find('.')-7:sampleFile.find('.')]
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
