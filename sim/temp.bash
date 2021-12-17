
btb_seq=~/repos/btb-seq
genomes=~/temp/bfast-genomes-3/
reads=~/temp/bfast-reads-3/
output=~/temp/sequence-4/
sudo rm -r $output
sudo python validator.py sequence $btb_seq $genomes $reads $output --light