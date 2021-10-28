import os
import subprocess
import glob
import json
import shutil

from compare_snps import analyse

reference_path = './Mycobacterium_bovis_AF212297_LT78304.fa'

def run(cmd, cwd=None):
    """ Run a command and assert it exits with a non-zero exit code
    """
    # TODO: store stdout to a file
    returncode = subprocess.run(cmd, cwd=cwd).returncode

    if returncode:        
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % cmd, returncode)


def simulate_genome(output_prefix, num_snps=16000):
    # TODO: make sure output path does not exist

    output_path = output_prefix+"simulated-genome/"
    os.makedirs(output_path, exist_ok=False)

    # TODO: store stdout to a file
    run([
        "simuG.pl",
        "-refseq", reference_path,
        "-snp_count", str(num_snps),
        "-prefix", output_path+"simulated"
    ])

def simulate_reads(
    genome_fasta, 
    base_directory,
    read_length=150,
    seed=1
):
    output_path = base_directory + 'simulated_reads/'
    os.makedirs(output_path, exist_ok=False)

    read_1 = output_path+"simulated_S1_R1_X.fastq"
    read_2 = output_path+"simulated_S1_R2_X.fastq"

    # TODO: store stdout to a file
    run([
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
    ])

    # TODO: use python library instead?
    run(["gzip", read_1, read_2])


def btb_seq(pipeline_directory, reads_directory, base_directory):
    results_directory = base_directory + 'pipeline/'
    os.makedirs(results_directory, exist_ok=True)

    run(["bash", "./btb-seq", reads_directory, results_directory], cwd=pipeline_directory)

def performance_test_btb_seq(btb_seq_path, results_path):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to where btb-seq code is stored
            results_path (str): Path to where performance test results are stored

        Returns:
            None
    """

    # Validate Input
    if os.path.isdir(results_path):
        raise Exception("Output results path already exists")
    os.makedirs(results_path)

    if not os.path.isdir(btb_seq_path):
        raise Exception("Pipeline code repository not found")  

    # Set results directories
    fasta_path = results_path + 'simulated-genome/simulated.simseq.genome.fa'
    reads_path = results_path + 'simulated_reads/'
    simulated_snps = results_path + "simulated-genome/simulated.refseq2simseq.map.txt"
    pipeline_snps = btb_seq_path + 'snpTables/simulated.tab'

    # Copy over the repo
    repo_backup_path=results_path+"btb-seq/"
    run(["cp", "-r", btb_seq_path, repo_backup_path])

    # Run Simulation
    simulate_genome(results_path)
    simulate_reads(fasta_path, results_path)    
    btb_seq(btb_seq_path, reads_path, results_path)

    # Analyse Results
    # HACK: this could easily break if additioanl files are present
    pipeline_directory = glob.glob(results_path + 'pipeline/*')[0] + '/'
    
    stats = analyse(simulated_snps, pipeline_snps)

    # Write output
    with open(results_path+"stats.json", "w") as file:
        file.write(json.dumps(stats, indent=4))

def main():
    pipeline_path = '/home/aaronfishman/repos/btb-seq/'
    results_path = '/home/aaronfishman/ebs/pipeline-results/btb-seq-6/'

    # TODO: argparse

if __name__ == '__main__':
    main()