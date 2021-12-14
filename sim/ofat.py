import os
import sys
import argparse
import glob
import json

import pandas as pd
import validator


"""
    Run the performance benchmarking tool against a multiple git branches
"""

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'

# Initial list of branches that we are testing
DEFAULT_BRANCHES = [
    "BaseQual",
    "MapQual",
    "Ploidy",
    "VarQual",
    "ReadDepth",
    "ForRevAltRead",
    "AltProportion",
    "SNPwindow",
    "RepeatMask",
    "v1",
    "v2"
]

def ofat(btb_seq_path, output_path, branches=DEFAULT_BRANCHES, 
         reference_path=DEFAULT_REFERENCE_PATH):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to btb-seq code is stored
            results_path (str): Output path to performance test results
            branches (list): List of strings for each git branch to test
    """
    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    output_path = os.path.join(output_path, '')     
    
    # generate samples
    samples = validator.standard_samples()
    # simulate reads
    validator.simulate(output_path, samples, reference_path)
    
    # Benchmark the branches
    for branch in branches:

        try:
            validator.performance_test(
                output_path,
                btb_seq_path,
                samples, 
                branch=branch,
                light_mode=True
            )
        except Exception as e:
            print(e)
            print(f"***FAILED BRANCH: {branch}****", branch)

    # Analyse results
    return None
    #return analyse(results_path)

def analyse(root_path):
    """ Analyse results from an ofat run

        Parameters:
            root_path (str): Path to parent directory where 
                results from ofat trials are stored
    """

    # Determine path of each trial
    paths = glob.glob(root_path + '/*')

    # Initialise output data columns
    fns = []
    fps = []
    tps = []
    branch_names = []

    # Collect data from each trial
    for path in paths:
        # Initialised trial data
        branch = os.path.basename(path)
        stats_path = path + '/stats.json'

        fn = 'FAIL'
        fp = 'FAIL'
        tp = 'FAIL'

        # Log 
        if os.path.exists(stats_path):
            with open(stats_path, 'r') as file:
                data = json.load(file)

            fn = data['fn']
            fp = data['fp']
            tp = data['tp']

        # Add to output
        fns.append(fn)
        fps.append(fp)
        tps.append(tp)
        branch_names.append(branch)

    return pd.DataFrame(data={
        "branch": branch_names,
        "fn": fns,
        "fp": fps,
        "tp": tps
    })

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Run the performance benchmarking tool against a number of git branches")

    parser.add_argument("btb_seq", help="path to btb-seq code")
    parser.add_argument("results", help="path to results directory")

    args = parser.parse_args(sys.argv[1:])

    # Run
    df = ofat(args.btb_seq, args.results)
    print("OFAT run completed, results:")
    print(df)