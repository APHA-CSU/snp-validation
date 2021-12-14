import subprocess
import pandas as pd
from io import StringIO

def run(cmd, *args, capture_output=False, **kwargs):
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
    ps = subprocess.run(cmd, *args, capture_output=capture_output, **kwargs)

    returncode = ps.returncode
    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % (cmd, returncode))

    if capture_output:
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
