import csv
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=True)
parser.add_argument('-o', '--output', help="the name of the output file", required=True)
parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population")
parser.add_argument('-s', '--synonymousPairing', action='store_true', help="include this flag to calculate the total number of synonymous codon pairs", required=False)
args = parser.parse_args()
allInputFiles = args.input
key = pd.read_csv(args.keyFile)

pops = {'GBR' : [], 'FIN' : [], 'CHS' : [], 'PUR' : [], 'CDX' : [], 'CLM' : [],
        'IBS' : [], 'PEL' : [], 'PJL' : [], 'KHV' : [], 'ACB' : [], 'GWD' : [],
        'ESN' : [], 'BEB' : [], 'MSL' : [], 'STU' : [], 'ITU' : [], 'CEU' : [],
        'YRI' : [], 'CHB' : [], 'JPT' : [], 'LWK' : [], 'ASW' : [], 'MXL' : [],
        'TSI' : [], 'GIH' : []}

codonDict = {'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 
          'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 
          'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 
          'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 
          'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 
          'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 
          'TTG': 0, 'TTT': 0}

aaDict = {'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 
          'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 
          'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, 'Z': 0}

if args.synonymousPairing:
    numColumns = 25
    myDict = aaDict
else:
    numColumns = 67
    myDict = codonDict
    
d = dict()
for x in allInputFiles:
    with open(x, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Sample':
                codons = row[3:]
            else:
                for i in range(3,numColumns):
                    myDict[codons[i - 3]] += int(row[i])
    sample = x[x.find('/') + 1 :-11]
    d[sample] = []
    for i in range(3,numColumns):
        d[sample].append(myDict[codons[i - 3]])
    myDict = myDict.fromkeys(myDict, 0)

key = pd.read_csv(args.keyFile)
header = ['Sample','Subpopulation','Superpopulation'] + codons
with open(args.output, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(header)
    for sample,totals in d.items():
        subpopulation = key[key['Sample'] == sample]['Subpopulation'].values[0]
        superpopulation = key[key['Sample'] == sample]['Superpopulation'].values[0]
        writer.writerow([sample] + [subpopulation] + [superpopulation] + totals)
