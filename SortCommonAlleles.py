import csv
import argparse

def formatRow(row):
    string_row = row[0] + ',' + row[1]
    if args.cotrna:
        for x in range(2,24):
            string_row += ',' + str(int(float(row[x])))
    else:
        for x in range(2,66):
            string_row += ',' + str(int(float(row[x])))
    return string_row

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="Input csv Files",nargs='*',action="store", dest="input", required=True)
    parser.add_argument("-o", "--output",help="Output File",action="store",dest="output", required=True)
    parser.add_argument("-c", "--cotRNA",help="flag to merge co-tRNA codon pairing CSV files", action="store_true",dest="cotRNA", required=False)
    args = parser.parse_args()
    return args

args = parseArgs()
allInputFiles = args.input
outFile = args.output

myDict = dict()
for x in allInputFiles:
    with open(x, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Sample':
                codons = row[3:]
            else:
                gene = formatRow(row[1:])
                if gene in myDict.keys():
                    myDict[gene].append(row[0])
                else:
                    myDict[gene] = list()
                    myDict[gene].append(row[0])

with open(outFile, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['Index', 'Gene', 'Name', 'Samples'] + codons)
    i = 1
    for key,value in myDict.items():
        gene_row = key.split(',')
        sampleList = str(value).replace("'",'').replace('[','').replace(']','')
        writer.writerow([i] + gene_row[0:1] + [gene_row[1]] + [sampleList] + gene_row[2:])
        i += 1