 blastn -task blastn -query families_nr.fasta   -subject IRGSP-1.0_genome.fasta -outfmt "6 qseqid sseqid qstart qend sstart send mismatch gaps pident evalue length qlen slen qcovs" > results.csv
