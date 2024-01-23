#########################################################################
# File Name: ldsc.sh
# 使用 annot 文件计算 LD 分数(以PT细胞marker基因为例）
#########################################################################
#!/usr/bin/bash
for j in {1..22}
do
  python 03.ldsc.envs/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${j} \
  --ld-wind-cm 1 --annot 02.make_annot.out/celltype_marker/PT.bed.${j}.annot.gz \
  --thin-annot  --out 04.ldsc.out/celltype_marker/PT.ldsc.${j}
done

