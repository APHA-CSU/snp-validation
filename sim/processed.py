class ProcessedSample:
    def __init__(self, simulated_genome, sequenced_sample):
        if simulated_genome.name != sequenced_sample.name:
            raise Exception("Sample names do not match")

        self.name = simulated_genome.name
        self.genome = simulated_genome
        self.sequenced = sequenced_sample

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
        raise Exception(f"""Genome sample names are different to sequenced sample names
            Genomes: {genome_dict.keys()}
            Sequenced: {sequenced_dict.keys()}
        """)
    
    # Match
    samples = []
    for name in genome_dict:
        sample = ProcessedSample(genome_dict[name], sequenced_dict[name])
        samples.append(sample)

    return samples
