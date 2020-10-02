
#add path/to/vcf.gz
VCFNAME=ancientA_filtered.vcf.gz
#add path/to/genome.fna
GENOME=~/Genomes/Bacillus_anthracis/ames_ancestor.fa

#index vcf and overwrite if indexed already
tabix -f $VCFNAME

#for each sample in the vcf (bcftools query output) use vcf-consensus to make a fasta file of each sample
for f in `bcftools query -l $VCFNAME`
do
bcftools consensus -f $GENOME -s $f -i 'type="snp"' $VCFNAME | sed 's/NC_007530_Bacillus_anthracis_Ames_Ancestor/'$f'/g' > $f".fna"
done

rm -rf renamed/
mkdir renamed

for f in *.fna; do awk '/^>/ {gsub(/.fna(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' $f > renamed/$f; done

rm *.fna
