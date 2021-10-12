from Bio import SeqIO

actual_filepath = '/home/aaronfishman/simulated-results/Results_ref_15Sep21/consensus/ref.fas'
reference_filepath = '/home/aaronfishman/repos/btb-seq/references/Mycbovis-2122-97_LT708304.fas'


def load_fasta(filepath):
    with open(filepath) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record._seq)

def hamming_distance(a, b):
    if len(a) != len(b):
        raise Exception('sequences not the same length %i and %i'%(len(a), len(b)))

    distance = 0

    for i, _ in enumerate(a):
        if a[i] != b[i]:
            distance += 1
    
    return distance


reference = load_fasta(reference_filepath)
actual = load_fasta(actual_filepath)
hamming = hamming_distance(reference, actual)


print("actual len", len(actual))
print("reference len", len(reference))
print("hamming", hamming)
