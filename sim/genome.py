import os
import glob

import utils

class SimulatedGenome:
    """ Simulated genome data generated by simuG  """

    def __init__(self, name, genome_path, snp_table_path, snp_vcf_path, indel_vcf_path):
        # Validate
        utils.assert_path_exists(genome_path)
        utils.assert_path_exists(snp_table_path)
        utils.assert_path_exists(snp_vcf_path)
        utils.assert_path_exists(indel_vcf_path)

        # Assign
        self.genome_path = genome_path
        self.snp_table_path = snp_table_path
        self.snp_vcf_path = snp_vcf_path
        self.indel_vcf_path = indel_vcf_path
        
        self.name = name

def from_directory(path):
    """ List of SimulatedGenome from directory """
    # Initialise
    postfix = '.simulated.simseq.genome.fa'
    path = os.path.join(path, '')
    filepaths = glob.glob(path + f'*{postfix}')

    # Construct
    samples = []
    for filepath in filepaths:
        name = os.path.basename(filepath[:-len(postfix)])
        
        genome_path = path + name + '.simulated.simseq.genome.fa'
        snp_table_path = path + name + '.simulated.refseq2simseq.map.txt'
        snp_vcf_path = path + name + '.simulated.refseq2simseq.SNP.vcf'
        indel_vcf_path = path + name + '.simulated.refseq2simseq.INDEL.vcf'

        samples.append(SimulatedGenome(name, genome_path, snp_table_path, snp_vcf_path, indel_vcf_path))

    return samples
