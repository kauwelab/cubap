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
    parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population")
    args = parser.parse_args()
    return args

def toCSV(rampDict, outFile):
    '''
    writes the dict to a columanr csv
    '''
    header = ['Gene', 'Name', 'Ramp Sequence', 'Superpopulation', 'Subpopulation', 'Sample']
    outFile = outFile + '.csv'
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for sample, geneList in rampDict.items():
            for gene in geneList:
                geneRow = [gene['Gene'], gene['Name'], gene['Ramp Sequence'], 
                           gene['Superpopulation'], gene['Subpopulation'], gene['Sample']]
                writer.writerow(geneRow)
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
        np = description[description.find("Name=") + 5 :
            description.find("Name=") + 5 +
            get_next_semicolon(description[description.find("Name=") + 5 : ])]
    except:
        np = ''
    return geneName, np

def parseGenes(sampleFile, sampleName, superpopulation, subpopulation):
    '''
    Sorts FastA file information into list
    '''
    genes = list()
    geneDict = {'Gene' : '', 'Name' : '', 'Ramp Sequence' : '', 'Superpopulation' : '', 'Subpopulation' : '', 'Sample' : ''}
    for seq_record in SeqIO.parse(sampleFile, "fasta"):
        seq = seq_record.seq
        description = seq_record.description
        if description.find('exception') == -1 and description.find('partial=true') == -1:
            geneDict['Gene'], geneDict['Name'] = get_gene_name(description)
            geneDict['Ramp Sequence'] = seq
            geneDict['Superpopulation'] = superpopulation
            geneDict['Subpopulation'] = subpopulation
            geneDict['Sample'] = sampleName
            genes.append(geneDict.copy())
    return genes

if __name__ =='__main__':
    '''
    Main.
    '''
    rampDict = dict()
    args = parseArgs()
    
    for sampleFile in args.input:
        sampleName = sampleFile[sampleFile.find('/') + 1: sampleFile.find('/') + 8]
        
        import pandas as pd
        df_pop = pd.read_csv(args.keyFile)
        superpopulation = df_pop[df_pop['Sample'] == sampleName.rstrip()]['Superpopulation'].item()
        subpopulation = df_pop[df_pop['Sample'] == sampleName.rstrip()]['Subpopulation'].item()
    
        rampDict[sampleName] = parseGenes(sampleFile, sampleName, superpopulation, subpopulation)
    
    if args.output is None:
        outFile = 'ramps_merged'
    else:
        outFile = args.output + 'ramps_merged'
        
    toCSV(rampDict, outFile)
