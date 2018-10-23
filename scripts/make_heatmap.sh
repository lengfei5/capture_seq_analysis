##########
### this script is to make heatmap and profile using BigWig files and bed files 
##########
cwd=`pwd`;
nb_cores=4;
OUT="$PWD/heatmaps_PartC/GC_correction"
DIR_peak="/groups/bell/jiwang/Projects/Jorge/PART_C/Peaks_Union_BEDs"

#Samples="H3K27me3 H2A119Ub"
Samples="Suz12 Ring1B Cbx7 Rybp"
Conditions="924E12 E12F01 515D10"
peak2use="10_10_PEAK_Union.bed"
binsize=200;

out0="WT_RybpKO_Eed_KO"
out_gcc=${out0}_peakExcluded_GCc

mkdir -p ${PWD}/logs

cd $OUT
for prot in $Samples
do 
    #prot="Ring1B";
    peak=`ls ${DIR_peak}/*"$peak2use" | grep $prot `;
    #ffs="`ls ${OUT}/*.bw | grep "$prot" |grep "AN312_" ` `ls ${OUT}/*.bw | grep "$prot" |grep "515D10_"` `ls ${OUT}/*.bw | grep "$prot" |grep "515D10H3_" `  `ls ${OUT}/*.bw | grep Input`"
    ffs=""
    for c in $Conditions
    do 
	#echo $c
	ffs="$ffs `ls ${OUT}/*.bw | grep "$prot" | grep $c `"
    done
    ffs="$ffs `ls ${OUT}/*.bw | grep Input`"

    echo $prot 
    echo $peak 
    #echo $ffs
    
    ffs1=`echo $ffs | tr ' ' '\n' |grep -v peakExcluded.bw `
    ffs2=" `echo $ffs | tr ' ' '\n' |grep peakExcluded.bw` `echo $ffs | tr ' ' '\n' | grep Input ` "
    echo $ffs1
    echo $ffs2
    
    matrix_name="${prot}_${out0}_matrix4heatmap.gz"
    matrix2save="${prot}_${out0}_matrix4save.txt"
    heatmap="${prot}_${out0}_heatmap_profile_RPKM"
    
    matrix_name2="${prot}_${out_gcc}_matrix4heatmap.gz"
    matrix2save2="${prot}_${out_gcc}_matrix4save.txt"
    heatmap2="${prot}_${out_gcc}_heatmap_profile_RPKM"
   
    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N heatmap "module load python/2.7.3; module load deeptools/2.5.0.1-python2.7.3; computeMatrix reference-point --referencePoint=center -S "$ffs1" -R $peak -a 5000 -b 5000 --outFileName "$matrix_name" --outFileNameMatrix "$matrix2save" --skipZeros --binSize=$binsize  --sortUsingSamples 1 --sortUsing=mean; plotHeatmap -m "$matrix_name" --colorMap Blues --sortRegions descend -o ${heatmap}.pdf; plotProfile -m "$matrix_name" -o ${heatmap}.png --plotHeight 12 --plotWidth 15 --colors blue green navy; " 
    
    #echo 
    qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N heatmap "module load python/2.7.3; module load deeptools/2.5.0.1-python2.7.3; computeMatrix reference-point --referencePoint=center -S "$ffs2" -R $peak -a 5000 -b 5000 --outFileName "$matrix_name2" --outFileNameMatrix "$matrix2save2" --skipZeros --binSize=$binsize  --sortUsingSamples 1 --sortUsing=mean; plotHeatmap -m "$matrix_name2" --colorMap Blues --sortRegions descend -o ${heatmap2}.pdf; plotProfile -m "$matrix_name2" -o ${heatmap2}.png --plotHeight 12 --plotWidth 15 --colors blue green navy; "
    
    #break;

done

cd $cwd 





 



