import os
import math
import scipy.stats as ss

input_dir = 'sliding_window'
output_dir = 'zscore.data.all'
num_lines = 10

def get_extremal_values(xs, n, direction):
    topN = []
    for i in range(len(xs)):
        topN.append( (xs[i], i) )
        if direction == 'highest':
            topN.sort(reverse=True)
        else:
            topN.sort()
        if len(topN) > n:
            topN = topN[:-1]
    return topN


def process_file(filename):
    infile = open(os.path.join(input_dir, filename), 'r')
    header = infile.readline().strip().split('\t')

    def get_index(colname):
        return header.index(colname)

    data = [line.strip().split('\t') for line in infile]
    for line in data:
        for i in range(1, len(line)):
            line[i] = float(line[i])

    data = [line for line in data if not any(map(math.isnan, line[1:]))]
    if len(data) == 0:
        return None

    def docolumn(colname):
        return [line[get_index(colname)] for line in data]

    cols = list(map(docolumn, header))
    datadict = dict(zip(header, cols))

    for key in header[1:]:
        datadict[key+'-z'] = list(ss.zscore(datadict[key], nan_policy="omit"))

    for key in header[1:]:
        pvals = [ss.norm.sf(x)*2 for x in datadict[key+'-z']]
        datadict[key+'-p'] = pvals


    datadict['low-high'] = [int(x.split('-')[0]) for x in datadict['low-high']]


    
    long_header = ['low-high', 'dnds', 'dn', 'ds', 'Pi', 'dnds-z', 'dn-z', 'ds-z', 'Pi-z', 'dnds-p', 'dn-p', 'ds-p', 'Pi-p']
        
    return datadict

def get_row(i, d):
    return [val[i] for key,val in d.items()]

def assemble_rows(dd):
    numrows = len(dd[list(dd.keys())[0]])
    rows = [get_row(i, dd) for i in range(numrows)]
    return rows

def tabify(row):
    return '\t'.join(map(str, row)) + '\n'

def write_output(dd, filename):
    if dd == None:
        return None
    print(filename)
    output_filename = os.path.join(output_dir, filename)
    outfile = open(output_filename, 'w')
    outfile.write(tabify(list(dd.keys())))
    for r in assemble_rows(dd):
        outfile.write(tabify(r))
    outfile.close()

    

def main():
    for file in os.listdir(input_dir):
        write_output(process_file(file), file)

        
            


    
