from glob import glob
import json
import argparse
import utils
import os

import config
import validator

def download_test_case(case, download_path):
    """ Downloads reads from AWS """

    with open(config.DEFAULT_JSON_PATH) as f:
        x =json.load(f)

    try:
        y = x[case]
    except KeyError:
        print(f"no case named {case}")
        
    for sample in y.keys():

        sample = y[sample]
        r1_uri = sample["r1_uri"]
        r2_uri = sample["r2_uri"]

        if os.path.exists(download_path + os.path.basename(r1_uri)):
            print("read already downloaded, skipping")
        else:
            cmd = f"aws s3 cp --dryrun {r1_uri} {download_path}"
            utils.run(cmd, shell=True)

        if os.path.exists(download_path + os.path.basename(r2_uri)):
            print("read already downloaded, skipping")
        else:
            cmd = f"aws s3 cp --dryrun {r2_uri} {download_path}"
            utils.run(cmd, shell=True)


def btb_seq(btb_seq_path, download_path, output_path):
    """ runs btb-seq on reads downloaded in download_test_case"""

    btb_out = output_path + "/btbseq"
    validator.sequence(btb_seq_path, download_path, btb_out)

def snps(output_path):
    """run snp-sites on consensus files, then runs snp-dists on the results"""

    # glob filenames
    megafasta = ""
    consensus_path = output_path + "/btbseq/Results*/consensus/"
    paths = glob(consensus_path + "*_consensus.fas")
    fasta_list = []
    for path in paths:
        with open(path) as file: fasta = file.read()
        fasta_list.append(fasta)
    megafasta = "\n".join(fasta_list)

    
    if os.path.exists(output_path + 'snps'):
        pass
    else:
        os.mkdir(output_path + 'snps')
        
    snpsout = output_path + 'snps/'
    with open(f'{snpsout}combined.fas', 'w') as file: 
        file.write(megafasta)

    
    # run snp sites 
    cmd = f'snp-sites {snpsout}combined.fas -o {output_path}snpsites.fas'
    utils.run(cmd, shell=True)
    # produce vcf of snps- for investigating sites
    cmd = f'snp-sites {snpsout}combined.fas -v -o {output_path}snpsites.vcf'
    utils.run(cmd, shell=True)

    # run snp-dists
    cmd = f'/home/cameronnicholls/snp-dists/snp-dists {output_path}snpsites.fas > {output_path}snps.tab'
    utils.run(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="download sample fastqs from aws from input of sample refs")
    
    parser.add_argument("case")
    parser.add_argument("download_path")
    parser.add_argument("output_path")
    parser.add_argument("btb_seq_path")
    args = parser.parse_args()
    
    download_test_case(args.case, args.download_path)
    btb_seq(args.btb_seq_path, args.download_path, args.output_path)
    snps(args.output_path)