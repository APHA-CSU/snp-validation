# Load a 

from sample import RandomSample, VcfSample
import validator

snp_vcf = '/home/aaronfishman/mnt/fsx-027/snippy/AF-21-07727-19.vcf'
results_path = '/home/aaronfishman/temp-results/fold-coverage-test-5/'
btb_seq = '/home/aaronfishman/repos/btb-seq'

genome_len = 4349904
read_len = 150

samples = [
    VcfSample(snp_vcf, num_read_pairs=0.5*10*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*11*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*12*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*13*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*14*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*15*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*16*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*17*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*18*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*19*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*20*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*30*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*40*genome_len/read_len),
    VcfSample(snp_vcf, num_read_pairs=0.5*50*genome_len/read_len),
]

validator.performance_test(results_path, btb_seq, samples=samples)