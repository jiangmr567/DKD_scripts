#########################################################################
# File Name: make_annot.sh
# 创建注释文件(以PT细胞marker基因为例）
#########################################################################
#!/usr/bin/bash
cd 01.celtype.DEG.star.end/celltype_marker
for j in {1..22}
do
  python 03.ldsc.envs/make_annot.py --bed-file PT.bed --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${j}.bim \
  --annot-file 02.make_annot.out/celltype_marker/PT.bed.${j}.annot.gz
done
