obtain MITEs
genome-wide search
filter results
extract FASTA
MSA same family
search for best alignments, check CN



md5sum Triticum_aestivum.IWGSC.dna.toplevel.fa 
773d7a2216e7e5dce0e27ee1e66d7780  Triticum_aestivum.IWGSC.dna.toplevel.fa

 blastn -evalue 10e-3 -query allmites.fasta -db db -outfmt "6 qseqid sseqid qstart qend sstart send mismatch gaps pident evalue length qlen slen qcovs" 