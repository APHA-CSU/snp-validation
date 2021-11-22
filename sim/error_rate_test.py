# Load a 

from sample import RandomSample, VcfSample
import validator

snp_vcf = '/home/aaronfishman/mnt/fsx-027/snippy/AF-21-07727-19.vcf'
results_path = '/home/aaronfishman/temp-results/per-base-error-test-3/'
btb_seq = '/home/aaronfishman/repos/btb-seq'

samples = [
    VcfSample(snp_vcf, per_base_error_rate="0"),
    VcfSample(snp_vcf, per_base_error_rate="0.001-0.01"),
    VcfSample(snp_vcf, per_base_error_rate="0.01-0.1"),
    VcfSample(snp_vcf, per_base_error_rate="0.1-1"),
]

validator.performance_test(results_path, btb_seq, samples=samples)