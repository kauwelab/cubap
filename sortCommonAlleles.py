import csv
import argparse

def formatRow(row):
    if args.cotrna:
        headerLength = 24
    elif args.nucleotides:
        headerLength = 6
    else:
        headerLength = 66
    string_row = row[0] + ',' + row[1]
        for x in range(2,headerLength):
            string_row += ',' + str(int(float(row[x])))
    return string_row

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",help="Input csv Files",nargs='*',action="store", dest="input", required=True)
    parser.add_argument("-o", "--output",help="Output File",action="store",dest="output", required=True)
    parser.add_argument("-c", "--cotrna",help="flag to merge co-tRNA codon pairing CSV files", action="store_true",dest="cotrna", required=False)
    parser.add_argument("-s", "--samples",help="include sample names in output CSV", action="store_true",dest="samples",required=False)
    parser.add_argument("-n", "--nucleotides",help="flag to merge nucleotide composition CSV files", action="store_true",dest="nucleotides",required=False)
    args = parser.parse_args()
    return args

args = parseArgs()
allInputFiles = args.input
outFile = args.output
print(args)
myDict = dict()
for x in allInputFiles:
    with open(x, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Sample':
                columns = row[3:]
            else:
                gene = formatRow(row[1:])
                if gene in myDict.keys():
                    myDict[gene].append(row[0])
                else:
                    myDict[gene] = list()
                    myDict[gene].append(row[0])

with open(outFile, 'w') as csv_file:
    writer = csv.writer(csv_file)
    if args.samples:
        writer.writerow(['Index', 'Gene', 'Name', 'Samples'] + columns)
        i = 1
        for key,value in myDict.items():
            gene_row = key.split(',')
            sampleList = str(value).replace("'",'').replace('[','').replace(']','')
            writer.writerow([i] + gene_row[0:1] + [gene_row[1]] + [sampleList] + gene_row[2:])
            i += 1
    else:
        writer.writerow(['Index', 'Gene', 'Name'] + columns)
        i = 1
        for key,value in myDict.items():
            gene_row = key.split(',')
            writer.writerow([i] + gene_row[0:1] + [gene_row[1]] + gene_row[2:])
            i += 1
