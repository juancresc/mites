genome=/home/juan/Desktop/juan/bio/data/IWGSC/42/Triticum_aestivum.IWGSC.dna.toplevel.fa
mites=data/rice/allmites.fasta.filtered.csv.gff
organism=wheat

python 2.extract_blast.py -g $genome -o $organism -m $mites

python 3.clusterize.py -g $genome -o $organism -m $mites

blastn -task blastn -query data/rice/2892.fa -db /home/juan/Desktop/juan/bio/data/genomes/rice/db -outfmt  "6 qseqid sseqid qstart qend sstart send mismatch gaps pident evalue length qlen slen qcovs" > data/rice/2892.genome.csv