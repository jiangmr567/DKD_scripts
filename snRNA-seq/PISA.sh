#########################################################################
# File Name: pisa.sh
# 对原始测序读数进行过滤和解复用
#########################################################################
a="AUTO_file"
b="soupx_result"
for file in `ls $a`
do
	mkdir -p $b"/"$file"/"ALL
	cd $a"/"$file
	cd $(ls)"/"02.cDNAAnno
    cat beads_stat.txt | cut -f1 > $b"/"$file"/"ALL"/"ALL.cell_stat.txt
	PISA count -@ 4 -tag CB -anno-tag GN -umi UR -outdir $b"/"$file"/"ALL -list $b"/"$file"/"ALL"/"ALL.cell_stat.txt final.bam
	cd $b"/"$file"/"ALL
	gunzip *
	mv features.tsv genes.tsv
	echo $file
done