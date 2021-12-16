import os
import glob

from vcf_sample import VcfSample
from random_sample import RandomSample

def vcf_samples(datasets, snippy_dir='/mnt/fsx-027/snippy/'):
    """ Returns a list of VcfSample objects built from data sets contained
        in the datasets parameter. Each object is instantiated with a 
        unique seed value

        Parameters:
            datasets (tuple (str)): elements describe the base directory name
                for each dataset
            snippy_dir (str): path of parent directory of datasets
    """
    vcf_filepaths = []
    for set_i in datasets:
        vcf_dir = snippy_dir + set_i
        if not os.path.isdir(vcf_dir):
            raise Exception("Predefined SNP directory not found")
            
        vcf_dir = os.path.join(vcf_dir, '')
        vcf_filepaths.extend(glob.glob(vcf_dir+'*.vcf'))
    
    # sort vcf_filepaths for using consistent seed values accross runs
    vcf_filepaths.sort()
    samples = []
    seed_value = 0
    for filepath in vcf_filepaths:
        seed_value+=1 # different seed value for each sample
        samples.append(VcfSample(filepath, seed=seed_value, per_base_error_rate="0.001-0.01"))    

    return samples

def quick_samples():
    """ Returns a list containing a single RandomSample object
    """
    return [RandomSample(seed=1)]

def random_samples():
    """ Returns a list containing two RandomSample objects 
    """
    return [RandomSample(seed=5, num_snps=0, num_indels=0, per_base_error_rate="0.001-0.01"),
            RandomSample(seed=1, per_base_error_rate="0.001-0.01"),
            RandomSample(seed=666, per_base_error_rate="0.001-0.01")]

def standard_samples():
    """ Returns a list of Sample objects for the standard dataset
    """
    samples = vcf_samples(('standard',)) + random_samples() 
    return samples

def aph_samples():
    """ Returns a list of Sample objects for the aph dataset
    """
    return vcf_samples(('aph',))

def zwyer_samples():
    """ Returns a list of Sample objects for the zwyer dataset 
    """
    return vcf_samples(('zwyer',))

def all_samples():
    """ Returns a list of Sample objects for all VCF samples and two
        random samples
    """
    samples = vcf_samples(('aph', 'zwyer')) + random_samples()
    return samples