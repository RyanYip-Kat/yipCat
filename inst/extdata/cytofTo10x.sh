matDir=$1
cd $matDir  && awk '{print $1"\t"$2"\tGene Expression"}' genes.tsv > features.tsv && rm genes.tsv && sed -i 's/\./\-/' barcodes.tsv && gzip *
