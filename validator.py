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
    simulate_genome(results_path)

    # read simulation -- chop up that genome

    # btb-seq

    # performance analysis

    return

if __name__ == '__main__':
    main()