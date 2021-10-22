import os
import subprocess


reference_path = './Mycobacterium_bovis_AF212297_LT78304.fa'
#
def simulate_genome(output_prefix, num_snps=16000):
    # TODO: make sure output path does not exist

    output_path = output_prefix+"simulated-genome/"
    os.makedirs(output_path, exist_ok=True)

    # TODO: store stdout to a file
    returncode = subprocess.run([
        "simuG.pl",
        "-refseq", reference_path,
        "-snp_count", str(num_snps),
        "-prefix", output_path+"simulated"
    ]).returncode

    if returncode:
        raise Exception("simuG failed with exit code ", returncode)

def simulate_reads(
    genome_fasta, 
    base_directory,
    read_length=150,
    seed=1
):
    output_path = base_directory + 'simulated_reads/'
    os.makedirs(output_path, exist_ok=True)

    read_1 = output_path+"simulated_S1_R1_X.fastq"
    read_2 = output_path+"simulated_S1_R2_X.fastq"

    # TODO: store stdout to a file
    returncode = subprocess.run([
        "wgsim",
        "-1", str(read_length),
        "-2", str(read_length),
        "-S", str(seed),
        "-r", "0",
        "-R", "0",
        "-X", "0",
        "-e", "0",
        genome_fasta,
        read_1,
        read_2
    ]).returncode

    if returncode:
        raise Exception("wgsim failed with exit code ", returncode)


def main():
    pipeline_path = '/home/aaronfishman/repos/btb-seq/'
    results_path = '/home/aaronfishman/validation-results/btb-seq/'

    # Housekeeping
    # make 
    # if os.path.isdir(results_path):
    #     raise Exception("Output results path already exists")
    # os.makedirs(results_path)

    if not os.path.isdir(pipeline_path):
        raise Exception("Pipeline code repository not found")  

    # simuG - simulate the genome
    # simulate_genome(results_path)

    # read simulation -- chop up that genome
    fasta_path = results_path + 'simulated-genome/simulated.simseq.genome.fa'
    simulate_reads(fasta_path, results_path)

    # btb-seq

    # performance analysis

    return

if __name__ == '__main__':
    main()