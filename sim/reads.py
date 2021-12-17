import os
import glob

class Reads:
    """ A pair of fastq reads """
    def __init__(self, read_1_path, read_2_path, name="unnamed"):
        if not os.path.exists(read_1_path):
            raise Exception("Could not find read ", read_1_path)

        if not os.path.exists(read_2_path):
            raise Exception("Could not find read ", read_2_path)

        self.read_1_path = read_1_path
        self.read_2_path = read_2_path
        self.name = name

def from_directory(path, read_1_postfix='_S1_R1_X.fastq.gz', read_2_postfix='_S1_R2_X.fastq.gz'):
    # Find reads
    path = os.path.join(path, '')

    read_1s = sorted(glob.glob(path + "*" + read_1_postfix))
    read_2s = sorted(glob.glob(path + "*" + read_2_postfix))

    if len(read_1s) != len(read_2s):
        raise Exception(f"Number of R1 reads different to R2 reads: ({len(read_1s)}, {len(read_2s)})")

    # Construct
    reads = []
    for read_1, read_2 in zip(read_1s, read_2s):
        name_1 = os.path.basename(read_1)[:-len(read_1_postfix)]
        name_2 = os.path.basename(read_2)[:-len(read_2_postfix)]

        if name_1 != name_2:
            raise Exception(f"Name mismatch, Could not pair reads: {name_1}, {name_2}")

        reads.append(Reads(read_1, read_2, name=name_1))

    return reads
