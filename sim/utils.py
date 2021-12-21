import subprocess
from io import StringIO
import os
import glob

import pandas as pd

def run(cmd, *args, **kwargs):
    """ Run a command and assert that the process exits with a non-zero exit code.
        See python's subprocess.run command for args/kwargs

        Parameters:
            cmd (list): List of strings defining the command, see (subprocess.run in python docs)
            cwd (str): Set surr

        Returns:
            None: if capture_output == False (default)
            process' stdout (str): if capture_output == True 
    """
    # TODO: store stdout to a file
    ps = subprocess.run(cmd, *args, **kwargs)

    returncode = ps.returncode
    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % (cmd, returncode))

    if "capture_output" in kwargs and kwargs["capture_output"]:
        return ps.stdout.decode().strip('\n')

def bcf_summary(filepath='/home/aaronfishman/temp/filtered.vcf'):
    # Query
    columns = ['POS', 'AD0', 'AD1', 'DP']
    text = subprocess.check_output(["bcftools", "query", "-f",
        '%POS, %INFO/AD{0}, %INFO/AD{1}, %INFO/DP\n',
        filepath
    ], text=True)
    
    # Construct data frame
    return pd.read_csv(StringIO(text), header=None, names=columns)

def checkout(repo_path, branch):
    run(["git", "checkout", str(branch)], cwd=repo_path)

def assert_path_exists(path):
    if not os.path.exists(path):
        raise Exception("Could not find path: ", path)

def names_consistent(x, y):
    """ True if the name properties between two lists are consistent and unique. False otherwise """

    x_dict = {X.name: X for X in x}
    y_dict = {Y.name: Y for Y in y}

    # Validate
    if len(x_dict) != len(x):
        return False

    if len(y_dict) != len(y):
        return False

    if set(x_dict.keys()) != set(y_dict.keys()):
        return False

    return True

def get_results_path(path):
    """ Return full results path name to match with btb-seq output """
    # TODO: handle when glob does not return a unique path
    return glob.glob(path + '/Results_*')[0] + '/'