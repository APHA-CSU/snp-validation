
READS_DIRECTORY=~/simulated-reads/
RESULTS_DIRECTORY=~/simulated-results/

# Read Simulation
REFERENCE=~/repos/btb-seq/references/Mycbovis-2122-97_LT708304.fas
READ_1=$READS_DIRECTORY/blah_S1_R1_X.fastq
READ_2=$READS_DIRECTORY/blah_S1_R2_X.fastq
READ_LEN=150
SEED=0

wgsim -1 $READ_LEN -2 $READ_LEN -S 0 -r 0 -R 0 -X 0 -e 0 $REFERENCE $READ_1 $READ_2
gzip $READS_DIRECTORY/*.fastq


# Run the pipeline
cd ~/repos/btb-seq/
./btb-seq $READS_DIRECTORY $RESULTS_DIRECTORY

