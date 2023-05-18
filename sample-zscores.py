import os
import numpy as np
import pandas as pd
import statistics
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter


input_dir = 'zscore.data.all'
output_dir = 'zscore.spectrum.output'



def read_file(filename):
    full_filename = os.path.join(input_dir, filename)
    df = pd.read_csv(full_filename, sep='\t', header=0)
    return df

def build_distribution(df):
    means = []
    for _ in range(1000):
        means.append( int(df.sample(10)['low-high'].mean())  )
    means.sort()
    return means

def get_high_loc(colname, df):
    return int(df.nlargest(10, colname)['low-high'].mean())

def get_low_loc(colname, df):
    return int(df.nsmallest(10, colname)['low-high'].mean())

def get_zscore(testvalue, distribution):
    mu = statistics.mean(distribution)
    sigma = statistics.stdev(distribution)
    z = (testvalue - mu) / sigma
    return z

def percentile(val, distribution):
    d = len(distribution)
    n = sum([1 for i in distribution if i <= val])
    return n/d


def process_file(filename):
    df = read_file(filename)
    if df.shape[0] < 20:
        return None
    dist = build_distribution(df)
    vars = ['dnds-z', 'dn-z', 'ds-z', 'Pi-z']
    high_locs = {k: get_high_loc(k, df) for k in vars}
    high_zscores = {k:get_zscore(v, dist) for k,v in high_locs.items()}
    high_pvals = {k:scipy.stats.norm.sf(abs(v))*2 for k,v in high_zscores.items()}

    # high_percentiles = {k:percentile(v, dist) for k,v in high_locs.items()}
    low_locs = {k: get_low_loc(k, df) for k in vars}
    low_zscores = {k:get_zscore(v, dist) for k,v in low_locs.items()}
    low_pvals = {k:scipy.stats.norm.sf(abs(v))*2 for k,v in low_zscores.items()}
    # low_percentiles = {k:percentile(v, dist) for k,v in low_locs.items()}

    of_name = os.path.join(output_dir, filename)
    outfile = open(of_name, 'w')
    outfile.write("field\tzscore\tpval\n")
    for item in vars:
        z = high_zscores[item]
        p = high_pvals[item]
        outfile.write(f'{item}\t{z}\t{p}\n')

    for item in vars:
        z = low_zscores[item]
        p = low_pvals[item]
        outfile.write(f'{item}\t{z}\t{p}\n')
    
    
    # for k, v in high_percentiles.items():
    #     outfile.write(f'{k}-high\t{v}\n')
    # for k, v in low_percentiles.items():
    #     outfile.write(f'{k}-low\t{v}\n')
                
def main():
    for filename in os.listdir(input_dir):
        print(filename)
        process_file(filename)

        
