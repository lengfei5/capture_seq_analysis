##########
# This script is to make heatmap and profile for ChIP-seq data 
# Inputs provided are : BigWg files and bed files for peaks or genes (or GTF file for genes) 
# In fact this script is too context-dependent, because for most of case, the heatmap and profiles
# are regrouped across conditions for the same chip
##########
nb_cores=4;
OUT="$PWD/heatmap_profiles_chipseq"
#DIR_peak="/groups/bell/jiwang/Projects/Jorge/PART_C/Peaks_Union_BEDs"
peak="/groups/bell/jiwang/annotations/mm10/Refseq/mm10_Refseq_Gene_downloadfromUCSC.bed"
DIR_bigwig="$PWD/bigWigs"
cwd=`pwd`;

#Samples="H3K27me3 H2A119Ub"
#Samples="Suz12 Ring1B Cbx7 Rybp"
#Conditions="924E12 E12F01 515D10"
#peak2use="10_10_PEAK_Union.bed"
binsize=200;

#out0="WT_RybpKO_Eed_KO"
#out_gcc=${out0}_peakExcluded_GCc

mkdir -p $OUT;
mkdir -p ${PWD}/logs

cd $OUT
for ffs in $DIR_bigwig/*.bw
do 
    #ffs="$ffs `ls ${OUT}/*.bw | grep Input`"
    #echo $prot 
    echo $peak 
    echo $ffs
    out=`basename "$ffs"`
    out=${out%.bw}
    echo $out
    #ffs1=`echo $ffs | tr ' ' '\n' |grep -v peakExcluded.bw `
    #ffs2=" `echo $ffs | tr ' ' '\n' |grep peakExcluded.bw` `echo $ffs | tr ' ' '\n' | grep Input ` "
    #echo $ffs1
    #echo $ffs2
    
    matrix_name="${out}_matrix4heatmap.gz"
    matrix2save="${out}_matrix4save.txt"
    heatmap="${out}_heatmap_profile_RPKM"
    
    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N heatmap "module load python/2.7.3; module load deeptools/2.5.0.1-python2.7.3; \
computeMatrix scale-regions -S "$ffs" -R $peak --regionBodyLength 10000 -b 10000 -a 2000 --outFileName "$matrix_name" --outFileNameMatrix "$matrix2save" --skipZeros --binSize=$binsize  --sortUsingSamples 1 --sortUsing=mean; plotProfile -m "$matrix_name" -o ${heatmap}.png --plotHeight 12 --plotWidth 15; " 
            
    #break;

done

cd $cwd 





 



