class ProcessedSample:
    def __init__(self, simulated_genome, sequenced_sample):
        self.genome = simulated_genome
        self.sample = sequenced_sample

        # TODO: other convenience functions?

def from_list(genomes, sequenced):
    genome_dict = {g.name: g for g in genomes}
    sequenced_dict = {s.name: s for s in sequenced}

    # Validate
    if len(genome_dict) != len(genomes):
        raise Exception("Genomes with non-unique names found")

    if len(sequenced_dict) != len(sequenced):
        raise Exception("Sequenced samples with non-unique names found")

    if set(genome_dict.keys()) != set(sequenced_dict.keys()):
        raise Exception("Genome sample names are different to sequenced sample names")
    
    # Match
    samples = []
    for name in genome_dict:
        sample = ProcessedSample(genome_dict[name], sequenced_dict[name])
        samples.append(sample)

    return samples