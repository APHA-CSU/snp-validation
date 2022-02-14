from glob import glob
import subprocess
import argparse
import os

def snps(input_path, output_path):
    # glob filenames
    megafasta = ""
    paths = glob(input_path + "*_consensus.fas")
    for path in paths:
        file = open(path)
        fasta = file.read()
        megafasta = megafasta + fasta + '\n'

    if os.path.exists(output_path + 'output'):
        os.chdir(output_path + 'output')
    else:
        os.mkdir(output_path + 'output')
        os.chdir(output_path + 'output')

    file = open(r'combined.fas', 'w')
    file.write(megafasta)
    file.close()

    # run snp sites 
    cmd = 'snp-sites combined.fas -o snpsites.fas'
    subprocess.run(cmd, shell=True)

    # produce vcf of snps- for investigating sites
    cmd = 'snp-sites combined.fas -v -o snpsites.vcf'
    subprocess.run(cmd, shell=True) 

    # run snp-dists
    cmd = 'snp-dists snpsites.fas > snps.tab'
    subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="download sample fastqs from aws from input of sample refs")
    
    parser.add_argument("input_path", default="./")
    parser.add_argument("output_path", default = "./")
    args = parser.parse_args()
    
    snps(args.input_path, args.output_path)