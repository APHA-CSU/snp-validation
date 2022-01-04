import glob
import os

import utils

class SequencedSample:
    """ output data generated by the btb-seq pipeline"""

    def __init__(self, name, vcf_path, filtered_bcf_path, snp_table_path, genome_path):
        # Validate
        utils.assert_path_exists(vcf_path)
        utils.assert_path_exists(filtered_bcf_path)
        utils.assert_path_exists(snp_table_path)
        utils.assert_path_exists(genome_path)

        # Assign
        self.vcf_path = vcf_path
        self.filtered_bcf_path = filtered_bcf_path
        self.snp_table_path = snp_table_path
        self.genome_path = genome_path
        self.name = name

def from_results_dir(results_dir):
    """ A list of samples conststructed from the btb-seq results directory """
    #Validate input
    utils.assert_path_exists(results_dir)

    # Sample names
    # TODO: Take from csv rather than consensus directory?
    paths = glob.glob(results_dir + '/consensus/*_consensus.fas')
    names = [os.path.basename(p)[:-14] for p in paths]

    # Construct samples
    samples = []
    for name in names:
        sample = SequencedSample(name, 
            results_dir + f'/vcf/{name}.vcf.gz',
            results_dir + f'/filteredBcf/{name}_filtered.bcf',
            results_dir + f'/snpTables/{name}_snps.tab',
            results_dir + f'/consensus/{name}_consensus.fas'
        )

        samples.append(sample)

    return samples
