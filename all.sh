genome=/home/juan/Desktop/juan/bio/data/IWGSC/42/Triticum_aestivum.IWGSC.dna.toplevel.fa

genome=/home/juan/Desktop/juan/bio/data/genomes/rice/IRGSP-1.0_genome.fasta
mites=data/rice/allmites.fasta.csv
data=data/rice

families=data/rice/families/
clusters=data/rice/clusters/
clusters_file=data/rice/clusters.csv
fasta_clusters=data/rice/clusters.fasta
fs_clusters_dir=data/rice/fs_cluster/
organism=rice

python blast_filter.py -i $mites -o ${mites}.filtered
python blast2gff.py -i ${mites}.filtered -n MITE -o ${mites}.filtered.gff3
#extract sequences
python gff2fastas.py -i ${mites}.filtered.gff3 -g $genome -o $families
python clusterize.py -f $families -c $clusters --csv $clusters_file --fasta $fasta_clusters -g $genome --data $data

python cluster_fs_compare.py --csv $clusters_file -f $families -c $fs_clusters_dir