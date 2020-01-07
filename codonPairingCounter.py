#! /usr/bin/env python
import sys
import csv
import argparse
import re
from Bio import SeqIO

def makeAllPossibleCodons(co_trna):
	'''
	Input: rna is a flag to specify if the sequence is DNA or RNA. co_trna is a flag to
		specify co_tRNA codon pairing (i.e., same amino acid formed)
	Returns a set of all 64 possible codons (DNA or RNA) or all 20 amino acids.
	'''
	if co_trna:
		return set(['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
	from itertools import product
	codons = product("ACGT",repeat=3)
	codonsComb = set()
	for c in codons:
		codonsComb.add("".join(c))
	return codonsComb

def parseArgs():
    '''
	Argument parsing is done.
	Required to have an input file.
	'''
    parser = argparse.ArgumentParser(description='Find Identical and co-tRNA codon pairing.')
    parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
    parser.add_argument('-o', '--output', type=str, help='the directory to which output files will be saved')
    parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
    parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
    parser.add_argument("-b",help="Both Identical and Co-tRNA codon pairing",action="store_true",dest="both", required=False)
    parser.add_argument('-p', '--population', help='the population of the current sample being processed',dest="population")
    parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population",dest="keyfile")
    parser.add_argument("-l",type=str, help="Codon Table. Default: Standard",action="store",dest="codon_table", default="Standard", required=False)
    parser.add_argument('--columnar', action='store_true', help='a flag to output results in a columnar csv', required=False)
    args = parser.parse_args()
    if not args.input and not args.inputDir:
        sys.stdout.write("You must supply an input file with either -i or -id\n")
        sys.exit()
    if args.co_trna and args.both:
        sys.stdout.write("You cannot use both the co_trna (-c) and both (-b) flags.\n")
        sys.exit()
    return args

def getPairs(seq,orderedCodons):
    '''
    Counts the codon pairs in a single gene sequence.
    Returns a dictionary of codon/amino acid keys and count values
    '''
    footprint = args.footprint
    pairCount = dict()
    codons = []
    dna_codons = []
    if args.both or args.co_trna:
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        sequence = Seq(seq,generic_dna)
        aa = str(sequence.translate(table=args.codon_table))
        codons = re.findall(".",aa)
    else:
        codons = re.findall("...",seq)
    if args.co_trna: #To ensure that identical codon pairing does not form the amino acid
        dna_codons = re.findall("...",seq)
    lastFound = dict() #key= codon, value= position of last found codon with pairing #For co-trna: key = codon (amino acid), value= dict() where key=dna_codon (codon) and value = last position of it
    for x in range(len(codons)):
        curCodon = codons[x]
        if not curCodon in orderedCodons:
            continue
        if not args.co_trna:
            if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint): #Must be >= because if footprint is 2 and AAA is found at positions 3 and 4, 4-3 =1, which is 1 less than footprint size.
                lastFound[curCodon] =x
                continue
        else:
            if not curCodon in lastFound:
                lastFound[curCodon] = dict()
                lastFound[curCodon][dna_codons[x]] =x
                continue
            closestPos = -100
            for key,value in lastFound[curCodon].items():
                if key == dna_codons[x]:
                    continue
                if value >closestPos:
                    closestPos = value
            if (x - closestPos) >= footprint:
                lastFound[curCodon][dna_codons[x]] =x
                continue
        if curCodon in pairCount.keys():
            pairCount[curCodon] += 1
        else:
            pairCount[curCodon] = 1
        if args.co_trna:
            lastFound[curCodon][dna_codons[x]] =x
            continue
        lastFound[curCodon] = x
    return pairCount

def readOneFileBioPython(inputFile):
    '''
    Reads one sample fasta file.
    Creates a dictionary of genes and their pairing.
    '''
    genes = dict()
    genes_A = dict()
    orderedCodons = makeAllPossibleCodons((args.co_trna | args.both))
    for seq_record in SeqIO.parse(inputFile, "fasta"):
        description = seq_record.description
        sequence = str(seq_record.seq)
        if description.find('exception') == -1 and description.find('partial=true') == -1 and len(sequence)%3==0:
            gene = description.split("gene=")[1].split(";")[0]
            name = description.split("Name=")[1].split(";")[0]
            geneName = gene + ' ' + name

            if geneName in genes.keys():
                genes_A[geneName] = getPairs(sequence,orderedCodons)
            else:
                genes[geneName] = getPairs(sequence,orderedCodons)
    return genes, genes_A

def toCSV(sampleName, population, genes, genes_A, outFile):
    '''
    writes the dict to a csv
    '''
    if args.co_trna:
        aminoAcids = ['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        if args.both:
            header = ['Sample', 'Population', 'Gene', 'Name'] + aminoAcids
            outFile = outFile + '_both.csv'
        elif not args.both:
            header = ['Sample', 'Population', 'Gene', 'Name'] + aminoAcids
            outFile = outFile + '_aa.csv'
        codonRow = []
        with open(outFile, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(header)
            for gene, codonDict in genes.items():
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                for codon in aminoAcids:
                    if codon in genes[gene].keys():
                        codonRow.append(genes[gene][codon])
                    else:
                        codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
            for gene, codonDict in genes_A.items():
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                for codon in aminoAcids:
                    if codon in genes_A[gene].keys():
                        codonRow.append(genes_A[gene][codon])
                    else:
                        codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
    else:
        codons = ['AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC',
                  'AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT',
                  'CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC',
                  'CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT',
                  'GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC',
                  'TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT',
                  'TTA','TTC','TTG','TTT']
        header = ['Sample', 'Population', 'Gene', 'Name'] + codons
        outFile = outFile + '_i.csv'
        codonRow = []
        with open(outFile, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(header)
            for gene, codonDict in genes.items():
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                for codon in codons:
                    if codon in genes[gene].keys():
                        codonRow.append(genes[gene][codon])
                    else:
                        codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
            for gene, codonDict in genes_A.items():
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                for codon in codons:
                    if codon in genes_A[gene].keys():
                        codonRow.append(genes_A[gene][codon])
                    else:
                        codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
    print(outFile,'created.')

def toColumnarCSV(sampleName, population, genes, genes_A, outFile):
    '''
    writes the dict to a columnar csv
    '''
    orderedCodons = makeAllPossibleCodons((args.co_trna | args.both))
    if args.co_trna and args.both:
        header = ['Sample', 'Population', 'Gene', 'Name', 'Amino Acid (including synonymous codons)', 'Count']
        outFile = outFile + '_both.csv'
    elif args.co_trna and not args.both:
        header = ['Sample', 'Population', 'Gene', 'Name', 'Amino Acid', 'Count']
        outFile = outFile + '_aa.csv'
    else:
        header = ['Sample', 'Population', 'Gene', 'Name', 'Codon', 'Count']
        outFile = outFile + '_i.csv'
    codonRow = []
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        # write genes
        for gene, codonDict in genes.items():
            for codon in orderedCodons:
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                if codon in genes[gene].keys():
                    codonRow.append(codon)
                    codonRow.append(genes[gene][codon])
                else:
                    codonRow.append(codon)
                    codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
        # write alleles
        for gene, codonDict in genes_A.items():
            for codon in orderedCodons:
                codonRow = [sampleName, population, gene.split(' ')[0], gene.split(' ')[1]]
                if codon in genes_A[gene].keys():
                    codonRow.append(codon)
                    codonRow.append(genes_A[gene][codon])
                else:
                    codonRow.append(codon)
                    codonRow.append(0)
                writer.writerow(codonRow)
                codonRow = []
    print(outFile,'created.')

if __name__ =='__main__':
    '''
    Main.
    '''
    args = parseArgs()
    if args.codon_table:
        try:
            if args.both or args.co_trna:
                from Bio.Seq import Seq
                from Bio.Alphabet import generic_dna
                sequence = Seq("ATG",generic_dna)
                aa = str(sequence.translate(table=args.codon_table))
        except Exception:
            sys.stderr.write("Codon Table Does Not Exist.\n")
            sys.exit()
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
    input = args.input[0]
    genes, genes_A = readOneFileBioPython(input)
    if args.columnar is True:
        toColumnarCSV(sampleName, population, genes, genes_A, outFile)
    else:
        toCSV(sampleName, population, genes, genes_A, outFile)
