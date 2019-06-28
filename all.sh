genome=/home/juan/Desktop/juan/bio/data/IWGSC/42/Triticum_aestivum.IWGSC.dna.toplevel.fa

#genome=data/IRGSP-1.0_genome.fasta
mites=data/wheat/mites_consensus.csv
families=data/wheat/families/
onlyflanking=data/wheat/ofs/
flanking=data/wheat/fs/
flanking_clusters=data/wheat/fs_clusters/
clusters=data/wheat/clusters/
clusters_file=data/wheat/clusters.csv
fasta_clusters=data/wheat/clusters.fasta
fs_clusters_dir=data/wheat/fs_cluster/

mkdir $families
mkdir $onlyflanking
mkdir $flanking
mkdir $fs_clusters_dir
mkdir $flanking_clusters
mkdir $clusters

#python blast_filter.py -i $mites -o ${mites}.filtered
#python blast2gff.py -i ${mites}.filtered -n MITE -o ${mites}.filtered.gff3
#python gff2fastas.py -i ${mites}.filtered.gff3 -g $genome -o $families
python clusterize.py --pid 0.93 --min 40 -f $families --clusters $clusters --csv $clusters_file --fasta $fasta_clusters -g $genome
#extract FS of interesting clusters
python 1cluster_extract_fs.py --with_fs $flanking -f $families --fsdir $onlyflanking --clusters $clusters --csv $clusters_file --fasta $fasta_clusters -g $genome
python cluster_fs_compare.py  --onlyflanking $onlyflanking --flankingcluster $flanking_clusters > data/res1.csv