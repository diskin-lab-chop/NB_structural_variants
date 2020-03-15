echo --- creating output directory ---
export basedir=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/CNV_analysis/GISTIC_SNP_GROUPS/HR_MNA/
mkdir $basedir
echo --- running GISTIC ---
export segfile=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/CNV_analysis/GISTIC_SNP_GROUPS/SNP_singleTumor_422_HR_MNA.seg
#export markersfile=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/SCA_project/GISTIC/markertab_V13.txt
export refgenefile=/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC/refgenefiles/hg19.mat
/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/gp_gistic2_from_seg -v 30 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -twoside 1 -brlen 0.98 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme -js 2 -rx 0
#/mnt/isilon/maris_lab/target_nbl_ngs/Gonzalo/install/GISTIC_2_0_23/gp_gistic2_from_seg -v 30 -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -twoside 1 -brlen 0.98 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme -js 2 -rx 0

# echo 'sh SNP_singleTumor_422_HR_MNA.sh' | qsub -cwd -pe smp 2 -l h_vmem=12G,m_mem_free=11G -N HR_MNA

