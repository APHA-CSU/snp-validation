from utils import run
import os

import config

class Sample:
    """ Sample data required for generating genomes and reads """

    per_base_error_rate = config.DEFAULT_PER_BASE_ERROR_RATE

    @property
    def name(self):
        return f"unnamed-{type(self).__name__}-{str(id(self))}"

    def simulate_genome(self):
        raise NotImplementedError("Please implement this method")

    def _simulate_genome_base(self, reference_path, simulated_genome_path, params):
        cmd = ["simuG.pl",
            "-refseq", reference_path,
            "-prefix", simulated_genome_path + ".simulated"
        ]
        cmd.extend(params)

        run(cmd)
    
    def simulate_reads(
        self,
        simulated_genome_path,
        simulated_reads_path,
        read_length=150,
        outer_distance=330,
        random_dna_probability=0.01,
        rate_of_mutations=0,
        indel_mutation_fraction=0,
        indel_extension_probability=0
    ):   
        """ Simulates raw illumina reads from genome fasta file 

            Parameters:
                simulated_genome_path (str): path to simulated genome fasta file
                simulated_reads_path (str): path to output simulated reads fastq files
                read_length (int): number of bases for each read
                seed (int): the seed value for simulation
                outer_distance (int):
                random_dna_probability (int):
                rate_of_mutations (int):
                indel_mutation_fraction (int):
                indel extension probability (int):

            Returns:
                None  
        """
        genome_fasta_path = simulated_genome_path

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
            "-o", "1",
            "-z", str(self.seed),
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

        return read_1, read_2