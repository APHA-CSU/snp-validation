
REFERENCE=~/repos/btb-seq/references/Mycbovis-2122-97_LT708304.fas
READ_1=reference_1.fastq
READ_2=reference_2.fastq

wgsim -1 300 -2 300 -r 0 -R 0 -X 0 -e 0 $REFERENCE $READ_1 $READ_2