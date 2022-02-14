from glob import glob
import subprocess
import argparse
import utils
import os

def snps(input_path, output_path):
    # glob filenames
    megafasta = ""
    paths = glob(input_path + "*_consensus.fas")
    fasta_list = []
    for path in paths:
        with open(path) as file: fasta = file.read()
        fasta_list.append(fasta)
    megafasta = "\n".join(fasta_list)

    
    if os.path.exists(output_path + 'output'):
        pass
    else:
        os.mkdir(output_path + 'output')
        
    snpsout = output_path + 'output/'
    with open(f'{snpsout}combined.fas', 'w') as file: 
        file.write(megafasta)

    
    # run snp sites 
    cmd = f'snp-sites {snpsout}combined.fas -o {snpsout}snpsites.fas'
    utils.run(cmd, shell=True)
    # produce vcf of snps- for investigating sites
    cmd = f'snp-sites {snpsout}combined.fas -v -o {snpsout}snpsites.vcf'
    utils.run(cmd, shell=True)

    # run snp-dists
    cmd = f'snp-dists {snpsout}snpsites.fas > {snpsout}snps.tab'
    utils.run(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="download sample fastqs from aws from input of sample refs")
    
    parser.add_argument("input_path", default="./")
    parser.add_argument("output_path", default = "./")
    args = parser.parse_args()
    
    snps(args.input_path, args.output_path)