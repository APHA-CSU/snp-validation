from glob import glob
import json
import argparse
import sys

sys.path.append('.')
import config
import validator
import utils

def download_test_case(case, download_path):
    """ Downloads reads from AWS """

    with open(config.DEFAULT_JSON_PATH) as f:
        test_cases =json.load(f)

    if not case in test_cases:
        raise Exception(f"no case named {case}")
        
    for sample in test_cases[case]:

        r1_uri = sample["r1_uri"]
        r2_uri = sample["r2_uri"]

        cmd = f"aws s3 sync {r1_uri} {download_path}"
        utils.run(cmd, shell=True)
       
        cmd = f"aws s3 sync {r2_uri} {download_path}"
        utils.run(cmd, shell=True)


def btb_seq(btb_seq_path, download_path, output_path):
    """ runs btb-seq on reads downloaded in download_test_case"""

    btb_out = output_path + "/btbseq"
    validator.sequence(btb_seq_path, download_path, btb_out)

def snps(output_path):
    """run snp-sites on consensus files, then runs snp-dists on the results"""

    # glob filenames
    consensus_path = output_path + "/btbseq/Results*/consensus/"
    paths = glob(consensus_path + "*_consensus.fas")

    snps_out = output_path
    with open(f'{snps_out}combined.fas', 'w') as snp_file: 

        for path in paths:

            with open(path) as fasta_file: fasta = fasta_file.read()
            snp_file.write(fasta + '\n')
    
    # run snp sites 
    cmd = f'snp-sites {snps_out}combined.fas -o {output_path}snpsites.fas'
    utils.run(cmd, shell=True)

    # produce vcf of snps- for investigating sites
    cmd = f'snp-sites {snps_out}combined.fas -v -o {output_path}snpsites.vcf'
    utils.run(cmd, shell=True)

    # run snp-dists
    cmd = f'snp-dists {output_path}snpsites.fas > {output_path}snps.tab'
    utils.run(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="download sample fastqs from aws from input of sample refs")
    
    parser.add_argument("case")
    parser.add_argument("download_path")
    parser.add_argument("output_path")
    parser.add_argument("btb_seq_path")
    args = parser.parse_args()

    download_test_case(args.case, args.download_path)
    #btb_seq(args.btb_seq_path, args.download_path, args.output_path)
    snps(args.output_path)