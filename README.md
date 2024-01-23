# Graduation Project: Diabetic kidney gene regulation mechanism revealed by Macaca fascicularis Model
## 各部分代码介绍
__data文件夹__：存放的是猴子基因注释及基因组注释，人类h19基因tss位点；  
__GWAS文件夹__：数据准备及可视化，LDSC三部曲：生成注释文件、进行LDSC打分及LDSC回归；  
__Integrated analysis文件夹__：snRNA-seq和snATAC-seq联合分析：peak2geneLink分析，TF调节因子潜在调控靶基因的鉴定，GRN调控网络的构建；  
__scJoint文件夹__：标签转移方法，参考scJoint官方教程（https://github.com/SydneyBioX/scJoint/tree/main/tutorial）；  
__snATAC-seqt文件夹__：单细胞核表观组下游分析:创建ArchR对象，识别差异峰，motif偏差富集，轨迹分析等；  
__snRNA-seq文件夹__：单细胞核转录组下游分析：PISA预处理，SoupX去除环境RNA污染，去除双胞，聚类注释，差异分析，cellchat细胞通讯，monocle3拟时序分析等； 
__utils文件夹__：相关运行必需的依赖工具、函数或公共代码，参考https://github.com/GreenleafLab/scScalpChromatin/blob/main； 
__species_comparison.R__：人猴鼠RNA-seq整合比较分析。
