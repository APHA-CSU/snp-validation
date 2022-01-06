import utils

class ProcessedSample:
    """ simulated genomes combined with and btb-seq sequence data 
        Useful for making comparisons when post-processing results
    """

    def __init__(self, simulated_genome, sequenced_sample):
        if simulated_genome.name != sequenced_sample.name:
            raise Exception("Sample names do not match")

        self.name = simulated_genome.name
        self.genome = simulated_genome
        self.sequenced = sequenced_sample

        # TODO: other convenience functions?

def from_list(genomes, sequenced):
    """ Construct processed samples from lists of SimulatedGenomes and SequencedSamples 
        Throws an error if the sample names are not consistent.
    """

    # Validate
    if not utils.names_consistent(genomes, sequenced):
        raise Exception("Genome names inconsistent with sequenced names")

    # Match
    genome_dict = {g.name: g for g in genomes}
    sequenced_dict = {s.name: s for s in sequenced}

    samples = []
    for name in genome_dict:
        sample = ProcessedSample(genome_dict[name], sequenced_dict[name])
        samples.append(sample)

    return samples
