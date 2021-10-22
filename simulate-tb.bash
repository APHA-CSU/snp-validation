ref=~/repos/btb-seq/references/Mycbovis-2122-97_LT708304.fas
# ref=./Mycbovis-2122-97_LT708304.fas.gz
ref=Mycbovis-2122-97_LT708304.fas.gz
ref=Mycbovis-2122-97_LT708304.fas.gz

ref=/home/aaronfishman/eles-reference/Mycobacterium_bovis_AF212297_LT78304.fa

perl ./simuG.pl \
    -refseq $ref \
    -snp_count 16000 \
    -prefix ~/testing-simuG/AF2122