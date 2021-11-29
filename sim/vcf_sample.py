import math
import os
import errno

import sample
from utils import run

class VcfSample(sample.Sample):
    def __init__(self, 
        predef_snp_path, 
        seed=1, 
        per_base_error_rate=sample.PER_BASE_ERROR_RATE,
        num_read_pairs = sample.NUM_READ_PAIRS,
    ):
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
        self.num_read_pairs = math.ceil(num_read_pairs)

        # TODO: check predef_snp_path exists
        self.predef_snp_path = predef_snp_path

    @property
    def name(self):
        # extract filename
        #TODO: this will fail if the extension is not 3 charachters long
        #TODO: this pattern is getting somewhat unweildy.
        #   etiher use something more general or save metadata separately
        basename = os.path.basename(self.predef_snp_path)[:-4]
        return f"{type(self).__name__}-{basename}-seed{self.seed}-error{self.per_base_error_rate}-readpairs-{self.num_read_pairs}"

    def simulate_genome(self, reference_path, simulated_genome_path):
        """ Simulate a genome fasta using simuG
            Simulated a genome with random SNPs

            TODO: rename simulated_genome_path to simulated_genome_prefix

            Parameters:
                reference_path (str): Path to reference genome
                simulated_genome_path (str): Path to simlated genome
                predef_snps_path (str): Path to snippy generated VCF file.
                seed (int): Seed value for simulation

            Returns:
                None
        """
        if not os.path.isfile(self.predef_snp_path):
            raise(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.predef_snp_path))
        
        # TODO better solution for storing decomposed VCF
        tmp_decomposed_vcf_path = simulated_genome_path + '.decomposed.vcf'
        
        # Decompose complex SNPs
        decompose_complex_snps(self.predef_snp_path, tmp_decomposed_vcf_path)
        params = ["-snp_vcf", tmp_decomposed_vcf_path, 
                #"-indel_vcf", predef_snp_path # tells simuG to also simulate predefined indels.
                "-seed", str(self.seed)]
        
        sample.simug(reference_path, simulated_genome_path, params)
        
        # clean tmp decomposed VCF
        run(['rm', tmp_decomposed_vcf_path])

def decompose_complex_snps(predef_snp_path, output_file_path):
    run(['vt', 
         'decompose_blocksub',
         predef_snp_path,
         '-o', output_file_path])
