from utils import run
import errno
import math 
import os

class Sample:

    # Do we need this?
    @property
    def name(self):
        # TODO: this is ugly, make this human-readable
        return f"unnamed-{str(id(self))}"

    def simulate_genome(self):
        raise NotImplementedError()

    def _simulate_genome_base(self, reference_path, simulated_genome_path, params):
        cmd = ["simuG.pl",
            "-refseq", reference_path,
            "-prefix", simulated_genome_path + ".simulated"]
        cmd.extend(params)
        run(cmd)
    
    def simulate_reads(
        self,
        simulated_genome_path,
        simulated_reads_path,
        read_length=150,
        seed=1,
        outer_distance=330,
        random_dna_probability=0.01,
        rate_of_mutations=0,
        indel_mutation_fraction=0,
        indel_extension_probability=0
    ):   

        genome_fasta_path = simulated_genome_path + self.name + '.simulated.simseq.genome.fa'

        output_prefix = simulated_reads_path + self.name
        
        run([
            "dwgsim",
            "-e", str(self.per_base_error_rate),
            "-E", str(self.per_base_error_rate),
            "-i",
            "-d", str(outer_distance),
            "-N", str(self.num_read_pairs),
            "-1", str(read_length),
            "-2", str(read_length),
            "-r", str(rate_of_mutations),
            "-R", str(indel_mutation_fraction),
            "-X", str(indel_extension_probability),
            "-y", str(random_dna_probability),
            "-H",
            "-z", str(seed),
            genome_fasta_path,
            output_prefix
        ])

        # How dwgsim chooses to name it's output fastq files
        dwgsim_read_1 = output_prefix + ".bwa.read1.fastq.gz"
        dwgsim_read_2 = output_prefix + ".bwa.read2.fastq.gz"
        
        # Output fastq filenames compatible with btb-seq
        read_1 = output_prefix + "_S1_R1_X.fastq.gz"
        read_2 = output_prefix + "_S1_R2_X.fastq.gz"
        
        # Rename output fastq files
        os.rename(dwgsim_read_1, read_1)
        os.rename(dwgsim_read_2, read_2)


class VcfSample(Sample):
    def __init__(self, 
                 predef_snp_path, 
                 seed=1, 
                 per_base_error_rate="0",
                 num_read_pairs = 144997,
                ):
        
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
        self.num_read_pairs = math.ceil(num_read_pairs)
        
        if os.path.isfile(predef_snp_path):
            self.predef_snp_path = predef_snp_path
        else:
            raise Exception("Cant Find the pipeline's snps table")
        
        # extract filename
        basename = os.path.splitext(os.path.basename(self.predef_snp_path))
        #TODO: this pattern is getting somewhat unweildy.
        #   etiher use something more general or save metadata separately
        self._name = f"{type(self).__name__}-{basename}-seed{self.seed}-error{self.per_base_error_rate}-readpairs-{self.num_read_pairs}"

    @property
    def name(self):
        return self._name

    def _decompose_complex_snps(self, output_file_path):
        run(['vt', 
            'decompose_blocksub',
            self.predef_snp_path,
            '-o', output_file_path])

    def simulate_genome(self, reference_path, simulated_genome_path, seed=1):
        """ Simulated a genome with random SNPs

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
        self._decompose_complex_snps(tmp_decomposed_vcf_path)
        params = ["-snp_vcf", tmp_decomposed_vcf_path, 
                "-indel_vcf", tmp_decomposed_vcf_path,
                "-seed", str(seed)]
        self._simulate_genome_base(reference_path, simulated_genome_path, params)
        # clean tmp decomposed VCF
        run(['rm', tmp_decomposed_vcf_path])

class RandomSample(Sample):
    def __init__(self, 
                 num_snps=16000, 
                 num_indels=1600, 
                 seed=1, 
                 per_base_error_rate="0",
                 num_read_pairs = 144997
                 ):
        
        self.num_snps = num_snps
        self.num_indels = num_indels
        self.seed = seed
        self.per_base_error_rate = per_base_error_rate
        self.num_read_pairs = math.ceil(num_read_pairs)
        
        self._name = f"{type(self).__name__}-snps{self.num_snps}-indels{self.num_indels}-seed{self.seed}"
    
    @property
    def name(self):
        return self._name

    def simulate_genome(self, reference_path, simulated_genome_path):
        """ Simulated a genome with random SNPs

            TODO: rename simulated_genome_path to simulated_genome_prefix
            Parameters:
                reference_path (str): Path to reference genome
                simulated_genome_path (str): Path to simlated genome
                num_snps (int): Number of random SNPs
                seed (int): Seed value for simulation

            Returns:
                None
        """
        params = ["-snp_count", str(self.num_snps),
                    "-indel_count", str(self.num_indels),
                    "-seed", str(self.seed)]
        self._simulate_genome_base(reference_path, simulated_genome_path, params)