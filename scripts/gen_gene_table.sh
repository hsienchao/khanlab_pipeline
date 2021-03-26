gtf=$1
if [ -z $gtf ];then
	gtf=/data/Clinomics/Ref/khanlab/GTF/gencode.v32.annotation.gtf
fi
echo -e "gene_id\tgene_type\tgene_name\tlevel" 
awk '$3=="gene"' $gtf | cut -f9 | perl -lane '{($gid,$gt,$gn,$level)=$_=~/gene_id "(.*?)"; gene_type "(.*?)"; gene_name "(.*?)"; (level .?);/;print "$gid\t$gn\t$gt\t$level"}'
