#########################################################################
# File Name: ldsc.regressions.sh
# ldsc回归分析(以PT细胞marker基因,eGFR的sumstats为例）
#########################################################################
#!/usr/bin/bash
python 03.ldsc.envs/ldsc.py  --h2-cts /sumstats/MW_eGFR.sumstats.gz \
 --ref-ld-chr 05.ldsc.regressions/1000G_EUR_Phase3_baseline/baseline. --n-blocks 1000 \
 --out 05.ldsc.regressions/03.regressions.result/celltype_marker/eGFR.regressions.out \
 --ref-ld-chr-cts 04.ldsc.out/celltype_marker/PT.ldsc.  \
 --w-ld-chr 05.ldsc.regressions/weights_hm3_no_hla/weights. --frqfile-chr 05.ldsc.regressions/1000G_Phase3_frq/1000G.EUR.QC.
