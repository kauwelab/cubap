import csv
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population")
args = parser.parse_args()
allInputFiles = args.input
key = pd.read_csv(args.keyFile)

pops = {'GBR' : [], 'FIN' : [], 'CHS' : [], 'PUR' : [], 'CDX' : [], 'CLM' : [],
        'IBS' : [], 'PEL' : [], 'PJL' : [], 'KHV' : [], 'ACB' : [], 'GWD' : [],
        'ESN' : [], 'BEB' : [], 'MSL' : [], 'STU' : [], 'ITU' : [], 'CEU' : [],
        'YRI' : [], 'CHB' : [], 'JPT' : [], 'LWK' : [], 'ASW' : [], 'MXL' : [],
        'TSI' : [], 'GIH' : []}

aaDict = {'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 
          'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 
          'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, 'Z': 0}

d = dict()
for x in allInputFiles:
    with open(x, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Sample':
                aas = row[3:]
            else:
                for i in range(3,25):
                    aaDict[aas[i - 3]] += int(row[i])
    sample = x[x.find('/') + 1 :-11]
    d[sample] = []
    for i in range(3,25):
        d[sample].append(aaDict[aas[i - 3]])
    aaDict = aaDict.fromkeys(aaDict, 0)
    #pops[key[key['Sample'] == sample]['Subpopulation'].values[0]].append(total)

key = pd.read_csv(args.keyFile)
header = ['Sample','Subpopulation','Superpopulation'] + aas
with open('total_co_tRNA_codon_pairing_freqs.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(header)
    for sample,totals in d.items():
        subpopulation = key[key['Sample'] == sample]['Subpopulation'].values[0]
        superpopulation = key[key['Sample'] == sample]['Superpopulation'].values[0]
        writer.writerow([sample] + [subpopulation] + [superpopulation] + totals)
