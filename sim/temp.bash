# 'bcftools query -f '%INFO/AD{0}, %INFO/AD{1}, %INFO/DP\n' /home/aaronfishman/temp/filtered.vcf

bcftools query -f '%POS, %INFO/AD{0}, %INFO/AD{1}, %INFO/DP\n' \
'/home/aaronfishman/temp/sites1/btb-seq-results/Results_simulated-reads_08Dec21/vcf/RandomSample-snps16000-indels1600-seed1.vcf.gz'
