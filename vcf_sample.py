import math
import os
import errno

from utils import run

from sample import Sample
from genome import SimulatedGenome
import config

class VcfSample(Sample):
    """ A sample generated from vcf file using simuG """
    def __init__(self, 
        predef_snp_path, 
        reference_path=config.DEFAULT_REFERENCE_PATH,
        seed=1, 
        per_base_error_rate=config.DEFAULT_PER_BASE_ERROR_RATE,
        num_read_pairs=config.DEFAULT_NUM_READ_PAIRS,
    ):

        # TODO: Validate
        self.reference_path = reference_path
        
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
        self.num_read_pairs = math.ceil(num_read_pairs)
        
        if os.path.isfile(predef_snp_path):
            self.predef_snp_path = predef_snp_path
        else:
            raise Exception("Cant find the pre-defined snp path")
        
    @property
    def name(self):
        # extract filename
        basename = os.path.splitext(os.path.basename(self.predef_snp_path))[0]
        #TODO: this pattern is getting somewhat unweildy.
        #   etiher use something more general or save metadata separately
        return f"{type(self).__name__}-{basename}-seed{self.seed}-error{self.per_base_error_rate}-readpairs-{self.num_read_pairs}"

    def decompose_complex_snps(self, output_file_path):
        run(['vt', 
            'decompose_blocksub',
            self.predef_snp_path,
            '-o', output_file_path])

    def simulate_genome(self, simulated_genome_path):
        """ Simulate the genome

            TODO: rename simulated_genome_path to simulated_genome_prefix

            Parameters:
                reference_path (str): Path to reference genome
                seed (int): Seed value for simulation

            Returns:
                None
        """

        # TODO: put the files into a temp directory

        if not os.path.isfile(self.predef_snp_path):
            raise(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.predef_snp_path))
        tmp_decomposed_vcf_path = simulated_genome_path + '.decomposed.vcf'
        
        # Decompose complex SNPs
        self.decompose_complex_snps(tmp_decomposed_vcf_path)
        params = ["-snp_vcf", tmp_decomposed_vcf_path, 
                  "-indel_vcf", tmp_decomposed_vcf_path,
                  "-seed", str(self.seed)]
        self._simulate_genome_base(self.reference_path, simulated_genome_path, params)
        
        # clean tmp decomposed VCF
        run(['rm', tmp_decomposed_vcf_path])

        # TODO: DRY with RandomSample
        snp_vcf_path = simulated_genome_path + '.simulated.refseq2simseq.SNP.vcf'
        indel_vcf_path = simulated_genome_path + '.simulated.refseq2simseq.INDEL.vcf'
        snp_table_path = simulated_genome_path + '.simulated.refseq2simseq.map.txt'
        genome_path = simulated_genome_path + '.simulated.simseq.genome.fa'

        return SimulatedGenome(self.name, genome_path, snp_table_path, snp_vcf_path, indel_vcf_path)
