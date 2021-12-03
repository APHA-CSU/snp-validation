import subprocess

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
    ps = subprocess.run(cmd, *args, **kwargs)

    if ps.returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % (cmd, ps.returncode))

    return ps