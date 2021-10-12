
REFERENCE=~/repos/btb-seq/references/Mycbovis-2122-97_LT708304.fas
READ_1=reference_1.fastq
READ_2=reference_2.fastq
READ_LEN=150

wgsim -1 $READ_LEN -2 $READ_LEN -r 0 -R 0 -X 0 -e 0 $REFERENCE $READ_1 $READ_2
gzip *.fastq