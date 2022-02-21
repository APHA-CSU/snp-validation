from glob import glob
import json
import argparse
import utils
import os

def download_test_case(case, download_path):
    with open('./test.json') as f:
        x =json.load(f)

        try:
            y = x[case]
        except:
            print(f"no case named {case}")
        
        num_samples = len(y)

        while num_samples > 0:
            sample = y[str(num_samples)]
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
            num_samples = num_samples - 1


def btb_seq(download_path, output_path):
    btb_out = output_path + "/btbseq"
    cmd = f" /home/cameronnicholls/Repos/btb-seq/btb-seq {download_path} {btb_out}"
    utils.run(cmd, shell=True, cwd = "/home/cameronnicholls/Repos/btb-seq/")


def snps(output_path):
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
    
    parser.add_argument("case", default="./")
    parser.add_argument("download_path")
    parser.add_argument("output_path", default = "./")
    args = parser.parse_args()
    
    download_test_case(args.case, args.download_path)
    btb_seq(args.download_path, args.output_path)
    snps(args.output_path)