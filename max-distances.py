import os
import math
import scipy.stats as ss
import pandas as pd

loc_source = "pi_merge_output"
zscore_source = "zscores.reduced"
output_dir = "zscores_finalized"



def get_short_filename(filename):
    fn = filename.split('.')
    return fn[0] + '.txt'

testfn = 'Escherichia_coli.txt.window.dnds.high'

def get_gene_loc_parameters(filename):
    shortfn = get_short_filename(filename)
    shortfn = os.path.join(loc_source, shortfn)
    df = pd.read_csv(shortfn, sep='\t', header=0)
    return math.floor(df['loc'].mean()), df['loc'].max()


def get_mean_loc(filename):
    fullpath = os.path.join(zscore_source, filename)
    df = pd.read_csv(fullpath, sep='\t', header=0)
    return math.floor(df['low-high'].mean())


def do_run(string):
    fullpath = os.path.join(output_dir, string[1:]+'.txt')
    output_file = open(fullpath, 'w')
    output_file.write('species\tori\toverall_mean\tnarrowed_mean\tterm\n')

    for f in [f for f in os.listdir(zscore_source) if f.endswith(string)]:
        sp = get_short_filename(f).split('.')[0]
        overall_mean, term = get_gene_loc_parameters(f)
        narrowed_mean = get_mean_loc(f)
        terms = map(str, [sp, '0', overall_mean, narrowed_mean, term])
        output_file.write('\t'.join(terms))
        output_file.write('\n')

    output_file.close()
    return
             

def main():
    strings = ['.'+x+y for x in ['dnds', 'dn', 'ds', 'Pi']
                       for y in ['.high', '.low']]
    for s in strings:
        do_run(s)


main()
