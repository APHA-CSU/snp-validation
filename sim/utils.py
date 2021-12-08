import subprocess
import pandas as pd
from io import StringIO

def run(cmd, *args, **kwargs):
    """ Run a command and assert that the process exits with a non-zero exit code.
        See python's subprocess.run command for args/kwargs

        Parameters:
            cmd (list): List of strings defining the command, see (subprocess.run in python docs)
            cwd (str): Set surr

        Returns:
            None
    """
    # TODO: store stdout to a file
    returncode = subprocess.run(cmd, *args, **kwargs).returncode

    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % (cmd, returncode))

def bcf_summary(filepath='/home/aaronfishman/temp/filtered.vcf', exclude=None):
    # excluded sites
    if exclude:
        exclusion_criteria = 'POS=='+' || POS=='.join(str(x) for x in exclude)
    else:
        exclusion_criteria = '1<0'
    
    # query BCF/VCF
    columns = ['POS', 'AD0', 'AD1', 'DP']
    text = subprocess.check_output(["bcftools", "query", "-f",
        '%POS, %INFO/AD{0}, %INFO/AD{1}, %INFO/DP\n',
        '-e', exclusion_criteria,
        filepath
    ], text=True)
    
    # construct data frame
    return pd.read_csv(StringIO(text), header=None, names=columns)
