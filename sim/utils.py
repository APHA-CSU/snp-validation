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

def bcf_summary(filepath='/home/aaronfishman/temp/filtered.vcf'):
      
    text = subprocess.check_output(["bcftools", "query", "-f",
        '%POS, %INFO/AD{0}, %INFO/AD{1}, %INFO/DP\n',
        filepath
    ], text=True)
    
    return pd.read_csv(StringIO(text), header=None,
        names = ['POS', 'AD0', 'AD1', 'DP']
    )
