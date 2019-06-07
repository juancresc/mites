for f in data/wheat/clusters_similar/*;
do 
    echo '************************'
    echo $f
    grep ">" $f | wc -l
    blastn -task blastn -query $f -subject data/db/trep-db_complete_Rel-16.fasta -outfmt  "6 sseqid"  | head -n 3      
done