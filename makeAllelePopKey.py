import csv
import argparse
import pandas as pd

def getPopCounts(index, samples):
    pops = {'Index': index, 'GBR' : 0, 'FIN' : 0, 'CHS' : 0, 'PUR' : 0, 'CDX' : 0, 'CLM' : 0,
        'IBS' : 0, 'PEL' : 0, 'PJL' : 0, 'KHV' : 0, 'ACB' : 0, 'GWD' : 0,
        'ESN' : 0, 'BEB' : 0, 'MSL' : 0, 'STU' : 0, 'ITU' : 0, 'CEU' : 0, 
        'YRI' : 0, 'CHB' : 0, 'JPT' : 0, 'LWK' : 0, 'ASW' : 0, 'MXL' : 0, 
        'TSI' : 0, 'GIH' : 0}
    for s in samples:
        pops[key[key['Sample'] == s]['Subpopulation'].values[0]] += 1
    return pops

def parseArgs():
    '''
    Argument parsing is done.
    Required to have an input file.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",help="Input Sorted File",action="store", dest="input", required=True)
    parser.add_argument("-o",help="Output File Name",action="store",dest="output", required=True)
    parser.add_argument('-k', '--keyFile', default='sampleKeyPop.csv', help="a csv key for identifying the sample's population", required=False)
    args = parser.parse_args()
    return args

if __name__ =='__main__':
    '''
    Main.
    '''
    args = parseArgs()
    inFile = args.input
    outFile = args.output
    
    key = pd.read_csv(args.keyFile)
    popCounts = list()
    
    with open(inFile, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] != 'Index':
                index = row[0]
                samples = list(row[3].split(', '))
                popCounts.append(getPopCounts(index, samples))
    
    header = ['Index', 'GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL',
              'PJL', 'KHV', 'ACB', 'GWD', 'ESN', 'BEB', 'MSL', 'STU', 'ITU',
              'CEU', 'YRI', 'CHB', 'JPT', 'LWK', 'ASW', 'MXL', 'TSI', 'GIH']
    
    with open(outFile, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for i in popCounts:
            writer.writerow([i['Index']] + [i['GBR']] + [i['FIN']] + [i['CHS']] + [i['PUR']] +
                            [i['CDX']] + [i['CLM']] + [i['IBS']] + [i['PEL']] + [i['PJL']] +
                            [i['KHV']] + [i['ACB']] + [i['GWD']] + [i['ESN']] + [i['BEB']] +
                            [i['MSL']] + [i['STU']] + [i['ITU']] + [i['CEU']] + [i['YRI']] +
                            [i['CHB']] + [i['JPT']] + [i['LWK']] + [i['ASW']] + [i['MXL']] +
                            [i['TSI']] + [i['GIH']])
