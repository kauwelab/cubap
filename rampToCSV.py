import csv
import argparse
import pandas as pd
import statistics

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input FASTA Files",nargs='*',action="store", dest="input", required=True)
parser.add_argument("-r", "--ramps", help="Input FASTA Files containing the ramp sequences",nargs='*',action="store", dest="ramps", required=False)
parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population")
parser.add_argument('-c', '--chart', action='store_true', help="output CSv only contains Gene Harmonic Mean RSCU, Ramp Harmonic Mean RSCU, Gene Length, Ramp Length", required=False)
parser.add_argument('-o', '--output', help="directory to which output files are saved")
args = parser.parse_args()
allInputFiles = args.input
allRampFiles = args.ramps
chart = args.chart
outputDir = args.output
df_pop = pd.read_csv(args.keyFile)

def make_column_number_list(d):
    longest = 0
    for ramp in d.values():
        if len(ramp) > longest:
            longest = len(ramp)
    nums = list(range(1,longest+1))
    return nums

def get_harmonic_means(ramp, rampLen):
    geneRSCU = statistics.harmonic_mean(ramp)
    rampRSCU = statistics.harmonic_mean(ramp[:rampLen])
    return geneRSCU, rampRSCU

def get_ramp_length(sample, gene, isoform, isAllele):
    first = False
    rampExists = False
    for file in allRampFiles:
        if sample in file:
            with open(file, 'r') as fp:
                line = fp.readline()
                while line:
                    if gene in line and isoform in line:
                        rampExists = True
                        if first == False:
                            first = True
                            firstRamp = fp.readline()
                        elif first == True:
                            first = False
                            secondRamp = fp.readline()
                    line = fp.readline()
    if rampExists == False:
        return 0
    else:
        if isAllele:
            return len(secondRamp)
        else:
            return len(firstRamp)

def to_csv_chart(sample, subpopulation, superpopulation, d, d_A):
    header = ['Sample','Subpopulation','Superpopulation', 'Gene', 'Name', 'Gene Harmonic Mean RSCU','Ramp Harmonic Mean RSCU','Gene Length','Ramp Length']
    with open(outputDir + sample + '_ramp_chart.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for gene,ramp in d.items():
            rampLen = get_ramp_length(sample, gene.split(' ')[0], gene.split(' ')[1], False)
            if rampLen != 0:
                geneRSCU, rampRSCU = get_harmonic_means(ramp, rampLen)
                writer.writerow([sample] + [subpopulation] + [superpopulation] + [gene.split(' ')[0]] + [gene.split(' ')[1]] + [str(geneRSCU)] + [str(rampRSCU)] + [str(len(ramp))] + [str(rampLen)])
        i = 0
        for gene,ramp in d_A.items():
            rampLen = get_ramp_length(sample, gene.split(' ')[0], gene.split(' ')[1], True)
            if rampLen != 0:
                i += 1
                geneRSCU, rampRSCU = get_harmonic_means(ramp, rampLen)
                writer.writerow([sample] + [subpopulation] + [superpopulation] + [gene.split(' ')[0]] + [gene.split(' ')[1]] + [str(geneRSCU)] + [str(rampRSCU)] + [str(len(ramp))] + [str(rampLen)])
def to_csv(sample, subpopulation, superpopulation, d, d_A):
    numList = make_column_number_list(d)
    if len(make_column_number_list(d_A)) > len(make_column_number_list(d)):
        numList = make_column_number_list(d_A)
    header = ['Sample','Subpopulation','Superpopulation'] + numList
    with open(outputDir + sample + '_rampies.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for gene,ramp in d.items():
             writer.writerow([sample] + [subpopulation] + [superpopulation] + ramp)
        for gene,ramp in d_A.items():
             writer.writerow([sample] + [subpopulation] + [superpopulation] + ramp)

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

d = dict()
d_A = dict()
setty = set()
is_description = False
for x in allInputFiles:
    with open(x, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            ramp = []
            if is_description == False and str(row[0])[0] == '>':
                is_description = True
                gene = get_gene_name(str(row))
            elif is_description == True and str(row[0])[0] != '>':
                if chart:
                    for rscu in row:
                        ramp.append(float(rscu.strip()))
                else:
                    for rscu in row:
                        if float(rscu.strip()) == 1.0:
                            ramp.append(int(float(rscu.strip())))
                        else:
                            ramp.append(round(float(rscu.strip()),3))
                is_description = False
                if gene in d.keys():
                    d_A[gene] = ramp
                else:
                    d[gene] = ramp
        sample = x[x.find('/') + 1 :-22]
        subpopulation = df_pop[df_pop['Sample'] == sample.rstrip()]['Subpopulation'].item()
        superpopulation = df_pop[df_pop['Sample'] == sample.rstrip()]['Superpopulation'].item()
        if chart:
            to_csv_chart(sample, subpopulation, superpopulation, d, d_A)
        else:
            to_csv(sample, subpopulation, superpopulation, d, d_A)
        print(sample)
