from utils import run
import errno
import math 
import os

# Default  dwgsim settings
READ_LEN = 150
OUTER_DISTANCE = 330
RANDOM_DNA_PROBABILITY = 0.01
RATE_OF_MUTATIONS = 0
INDEL_MUTATION_FRACTION = 0,
INDEX_EXTENSION_PROBABILITY = 0,
PER_BASE_ERROR_RATE = "0.001-0.01"
NUM_READ_PAIRS = 144997

class Sample:
    per_base_error_rate = PER_BASE_ERROR_RATE

    def simulate_genome(self):
        raise NotImplementedError()

    @property
    def name(self):
        # TODO: this is ugly, make this human-readable
        return f"unnamed-{str(id(self))}"

    def simulate_reads(self, 
        simulated_genome_dir, 
        simulated_reads_path, 
        seed=1,
    ):
        fasta_path = simulated_genome_dir + self.name + '.simulated.simseq.genome.fa'

        # How dwgsim chooses to name it's output fastq files
        output_prefix = simulated_reads_path + self.name
        dwgsim_read_1 = output_prefix + ".bwa.read1.fastq.gz"
        dwgsim_read_2 = output_prefix + ".bwa.read2.fastq.gz"

        # Output fastq filenames compatible with btb-seq
        read_1 = output_prefix + "_S1_R1_X.fastq.gz"
        read_2 = output_prefix + "_S1_R2_X.fastq.gz"

        run([
            "dwgsim",
            "-e", str(self.per_base_error_rate),
            "-E", str(self.per_base_error_rate),
            "-i",
            "-d", str(OUTER_DISTANCE),
            "-N", str(self.num_read_pairs),
            "-1", str(READ_LEN),
            "-2", str(READ_LEN),
            "-r", str(RATE_OF_MUTATIONS),
            "-R", str(INDEL_MUTATION_FRACTION),
            "-X", str(INDEX_EXTENSION_PROBABILITY),
            "-y", str(RANDOM_DNA_PROBABILITY),
            "-H",
            "-z", str(seed),
            fasta_path,
            output_prefix
        ])

        # Rename output fastq files
        os.rename(dwgsim_read_1, read_1)
        os.rename(dwgsim_read_2, read_2)

def simug(reference_path, simulated_genome_path, params):
    cmd = ["simuG.pl",
        "-refseq", reference_path,
        "-prefix", simulated_genome_path + ".simulated"
    ]

    cmd.extend(params)
    
    run(cmd)