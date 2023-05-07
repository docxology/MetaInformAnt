GO annotation based on eggNOGmapper


#install eggNOGmapper via conda

conda install eggnog-mapper

#download the supporting database
#http://eggnog5.embl.de/download/emapperdb-5.0.2/

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz


wget http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz #optional
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz #optional


#run eggnog-mapper

emapper.py --cpu 60 -i {sequence file, e.g. bee.fas} --output {output file name} -d euk -m diamond

#if diamond is not installed, use conda install -c bioconda diamond