import glob
import os

class SequencedSample:
    def __init__(self, name, vcf_path, filtered_bcf_path, snp_table_path):
        # Validate
        self.assert_path_exists(vcf_path)
        self.assert_path_exists(filtered_bcf_path)
        self.assert_path_exists(snp_table_path)

        # Assign
        self.vcf_path = vcf_path
        self.filtered_bcf_path = filtered_bcf_path
        self.snp_table_path = snp_table_path
        self.name = name

    def assert_path_exists(self, path):
        if not os.path.exists(path):
            raise Exception("Could not find path: ", path)

def from_results_dir(results_dir):
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
            results_dir + f'/snpTables/{name}_snps.tab'
        )

        samples.append(sample)

    return samples