obtain MITEs
genome-wide search
filter results
extract FASTA
MSA same family
search for best alignments, check CN



md5sum Triticum_aestivum.IWGSC.dna.toplevel.fa 
773d7a2216e7e5dce0e27ee1e66d7780  Triticum_aestivum.IWGSC.dna.toplevel.fa

 blastn -evalue 10e-3 -query allmites.fasta -db db -outfmt "6 qseqid sseqid qstart qend sstart send mismatch gaps pident evalue length qlen slen qcovs" 




Genome-wide analyses of miniature inverted-repeat transposable elements reveals new insights into the evolution of the Triticum-Aegilops group
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204972

Genome-Wide Analysis of Stowaway-Like MITEs in Wheat Reveals High Sequence Conservation, Gene Association, and Genomic Diversification
https://doi.org/10.1104/pp.112.204404

Transposable elements generate population-specific insertional patterns and allelic variation in genes of wild emmer wheat (Triticum turgidum ssp. dicoccoides)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5659041/

Mobilization of Stowaway-like MITEs in newly formed allohexaploid wheat species
https://link.springer.com/article/10.1007/s11103-012-9957-3


./bin/vsearch-2.13.4/bin/vsearch \
--cluster_fast data/wheat/allmites.fasta \
--threads 1 \
--strand both \
--clusters data/wheat/mites/ \
--iddef 1 \
--id 0.85