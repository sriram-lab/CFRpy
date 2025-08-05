"""
Python script to run `cfr_copt` across multiple conditions in parallel.

Example usage in command line: 
    > python apply_cfr_parallel.py --model_url <url> --data_path <path/to/data> --sheet_name <sheet_name> --thresh -1 1

Example usage with kwargs: 
    > python apply_cfr_parallel.py ... pfba_flag=False solver=glpk
"""

# Dependencies
import sys
import pickle
# import argparse
import pandas as pd
from tqdm import tqdm
from os import remove
from functools import partial
from itertools import repeat
from cfr import cfr_copt, Result
from multiprocessing import Pool
from cobra.io import load_matlab_model
from urllib.request import urlretrieve

'''# Define inputs from command line
parser = argparse.ArgumentParser(description='Inputs for apply_cfr')
parser.add_argument('--model_url', action='store', dest='url', required=True)
parser.add_argument('--data_path', action='store', dest='data_path', required=True)
parser.add_argument('--sheet_name', action='store', dest='sheet_name', required=False, default=None)
parser.add_argument('--thresh', action='store', dest='thresh', required=False, default=[-2, 2], nargs=2, type=float)
cmd_args = parser.parse_args()
url, data_path, sheet_name, thresh = cmd_args.url, cmd_args.data_path, cmd_args.sheet_name, cmd_args.thresh'''

# Define processing methods
def run_apply_async_multiprocessing(func, argument_list):

    with Pool() as pool: 
        jobs = [pool.apply_async(func=func, args=(*argument,)) if isinstance(argument, tuple) else pool.apply_async(func=func, args=(argument,)) for argument in argument_list]
    result_list_tqdm = []
    for job in tqdm(jobs):
        result_list_tqdm.append(job.get())

    return result_list_tqdm

def run_imap_multiprocessing(func, argument_list):

    with Pool() as pool: 
        result_list_tqdm = []
        for result in tqdm(pool.imap(func=func, iterable=argument_list), total=len(argument_list)):
            result_list_tqdm.append(result)

    return result_list_tqdm

def main(): 
    # Define ID
    id = 'eco'
    thresh = (-1, 1)

    # Load omics data
    if id=='eco':  
        url = 'http://bigg.ucsd.edu/static/models/iJO1366.mat' 
    elif id=='mtb': 
        url = 'http://bigg.ucsd.edu/static/models/iEK1008.mat'
    elif id=='human': 
        url = 'https://github.com/SysBioChalmers/Human-GEM/raw/refs/heads/main/model/Human-GEM.mat'

    # Define function inputs
    urlretrieve(url, './tmp_gem.mat')
    model = load_matlab_model('./tmp_gem.mat')
    data = pd.read_excel('data/test_cfr.xlsx', sheet_name=id, index_col=0, engine='openpyxl')

    # Process data
    genes = [gene.id for gene in model.genes]
    data = data[data.index.isin(genes)]

    # Determine up- and down-regulated genes
    off_dict = {col: data[data[col] < thresh[0]].index.tolist() for col in data.columns}
    on_dict = {col: data[data[col] > thresh[1]].index.tolist() for col in data.columns}

    # Define function arguments
    args = zip(repeat(model), on_dict.values(), off_dict.values())
    kwargs = {x.split('=')[0]: x.split('=')[1] for x in sys.argv if '=' in x}

    # Run function in parallel
    solutions = run_imap_multiprocessing(func=partial(cfr_copt, **kwargs), argument_list=args)

    # Define outputs
    sol_dict = {col: solutions[i] for i, col in enumerate(data.columns)}
    results = {col: Result(on_list=on_dict[col], off_list=off_dict[col], solution=sol_dict[col]) for col in data.columns}
    fname = model.id + '_test_cfr_output.pkl'
    with open(fname, 'wb') as f: 
        pickle.dump(results, f)

    # Run `cfr_copt` in parallel
    if __name__ == '__main__':
        main()