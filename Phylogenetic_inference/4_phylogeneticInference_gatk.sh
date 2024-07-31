PREFIX=final
PHYL_DIR=/scratch/projects/nib00015/phylogenomicsProject/gatk/phylogeneticInference
VCF_DIR=/scratch/projects/nib00015/phylogenomicsProject/gatk
SCRIPTS_DIR=/home/nibtve93/scripts/phylogeneticInference

AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py
ASCBIAS=/home/nibtve93/software/raxml_ascbias/ascbias.py

mkdir -p $PHYL_DIR/logFiles
mkdir -p $PHYL_DIR/raxml
mkdir -p $PHYL_DIR/iqtree
mkdir -p $PHYL_DIR/svdq


##################################
######Prepare for alignments######
##################################
for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.3 0.35 0.4 0.1 0.15 0.2 0.6 0.65 0.7 #
do
sbatch --output=$PHYL_DIR/logFiles/prepareForRaxml_maxmiss$i.oe $SCRIPTS_DIR/prepareForRaxml_SLURM.sh $VCF_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.vcf $PHYL_DIR
done


##################################
######Running RAxML for SNPs######
##################################
#run without bootstraps
for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7 #0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.3 0.35 0.4
do
NC=70
OUTGROUP="Cheirogaleusmaj_106245,Cheirogaleusmed_106354,Cheirogaleuscro_106189"
#run fast raxml algorithm
mkdir -p $PHYL_DIR/raxml/maxmiss$i
STAMATAKIS1=$(awk '{print $1}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/raxml-phylo_maxmiss$i.oe $SCRIPTS_DIR/raxml-phylo.sh $NC $OUTGROUP $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/raxml/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4
done

#run properly with bootstraps
for i in 0.05 0.25 #0.5 0.75 0.95 
do
NC=70
OUTGROUP="Cheirogaleusmaj_106245,Cheirogaleusmed_106354,Cheirogaleuscro_106189"
#run fast raxml algorithm (that does 20 unconstrained searches and 30 fast bootstrap replicates)
mkdir -p $PHYL_DIR/raxml/maxmiss$i
STAMATAKIS1=$(awk '{print $1}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/raxml-ng_maxmiss$i.boot.oe $SCRIPTS_DIR/raxml-ng.sh $NC $OUTGROUP $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/raxml/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.boot.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4
done

#speedup test for 0.05 threshold
for i in 1 4 8 12 16 20 40 76 #80 not tested because of overdetection
do
NC=$i
OUTGROUP="Cheirogaleusmaj_106245,Cheirogaleusmed_106354,Cheirogaleuscro_106189"
#run fast raxml algorithm (that does 20 unconstrained searches and 30 fast bootstrap replicates)
mkdir -p $PHYL_DIR/raxml/maxmiss0.05_speedup
STAMATAKIS1=$(awk '{print $1}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.noinv.phy.stamatakis)
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/raxml-ng_maxmiss0.05_speedup$i.oe $SCRIPTS_DIR/raxml-phylo.sh $NC $OUTGROUP $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.noinv.phy $PHYL_DIR/raxml/maxmiss0.05_speedup/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.noinv.speedup$i.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4
done


###################################
######Running IQ-TREE for SNPs#####
###################################
for i in 0.05 0.25
do
NC=80
UFBOOT=1000
#run iqtree
mkdir -p $PHYL_DIR/iqtree/maxmiss$i
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/iqtree-phylo_maxmiss$i.oe $SCRIPTS_DIR/iqtree-phylo.sh $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/iqtree/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv $UFBOOT
done

#test: without -alrt and -bnni
for i in 0.05 
do
NC=80
UFBOOT=1000
#run iqtree
mkdir -p $PHYL_DIR/iqtree/maxmiss$i
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/iqtree-phylo_no_alrt_bnni_maxmiss$i.oe $SCRIPTS_DIR/iqtree-phylo_no_alrt_bnni.sh $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/iqtree/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.no.alrt.bnni $UFBOOT
done

#test: with only -alrt 1000
for i in 0.05 
do
NC=80
UFBOOT=1000
#run iqtree
mkdir -p $PHYL_DIR/iqtree/maxmiss$i
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/iqtree-phylo_alrt_maxmiss$i.oe $SCRIPTS_DIR/iqtree-phylo_alrt.sh $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/iqtree/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.alrt $UFBOOT
done

#test: with only -bnni (increases computation time by about x2; decreases bootstrap values)
for i in 0.05
do
NC=80
UFBOOT=1000
#run iqtree
mkdir -p $PHYL_DIR/iqtree/maxmiss$i
sbatch --qos=7d -c $NC --output=$PHYL_DIR/logFiles/iqtree-phylo_bnni_maxmiss$i.oe $SCRIPTS_DIR/iqtree-phylo_bnni.sh $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/iqtree/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.bnni $UFBOOT
done


#run with checkpointing for thresholds 0.5, 0.75 and 0.95
NJOBS=12 #number of checkpointing runs per independent bpp run
UFBOOT=1000

for j in 0.95 #0.5 0.75 0.95
do
OUT_DIR=$PHYL_DIR/iqtree/maxmiss$j
mkdir -p $OUT_DIR

for i in $(seq 1 $NJOBS)
do
#logging
echo "submitting job $i"
#special case: if it's the first script
if [ $i == 1 ]
then
#declare dir for saving checkpoint files to (it will be created in gphocs script)
CHECK_OUT=$OUT_DIR/checkpoint_$i
#submit job and save printed text to variable
SUBTEXT=$(sbatch --output=$OUT_DIR/$i.iqtree.oe --account=nib00015 $SCRIPTS_DIR/iqtree_checkpoint.sh $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$j.noinv.phy $PHYL_DIR/iqtree/maxmiss$j/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$j.noinv $UFBOOT $CHECK_OUT)
#extract sbatch ID and store in variable
declare RUNID_${i}=${SUBTEXT##* }

#otherwise: submit rerun scripts
else
#assign input dir (which is output dir of previous iteration)
CHECK_IN=$CHECK_OUT
#declare dir for saving checkpoint files to (it will be created in gphocs script)
CHECK_OUT=$OUT_DIR/checkpoint_$i
#get name of variable that stored sbatch ID in previous iteration of the loop
VARNAME=RUNID_$(( ${i} - 1 ))
#submit rerun and save again printed text to variable
echo ${!VARNAME}
SUBTEXT=$(sbatch --output=$OUT_DIR/$i.iqtree.oe --account=nib00015  --dependency=afterany:${!VARNAME} $SCRIPTS_DIR/iqtree_checkpoint_rerun.sh $CHECK_IN $CHECK_OUT)
#see above
declare RUNID_${i}=${SUBTEXT##* }
fi
done

done

--account=nib00015 

#######################################################
####final run with bootstraps and 20 starting trees####
####not done									   ####
#######################################################


#####final run with bootstraps and 20 starting trees
#####do on nibur64b!!!!!####################
SCRIPTS_DIR=/home/nibur64b/scripts/phylogenomicsProject
PHYL_DIR=/scratch/usr/nibur64b/phylogenomicsproject/phylogeneticInference
NC=75
OUTGROUP="Cheirogaleusmaj_106245,Cheirogaleusmed_106354,Cheirogaleuscro_106189"
mkdir -p $PHYL_DIR/logFiles
#1) for maxmiss 0.25 run 10 times with parsimony starting tree and 10 times with random starting tree (no checkpointing necessary here because of fast runs); for 0.05 not necessary because already done on nibtve93 server
i=0.25
mkdir $PHYL_DIR/maxmiss$i
STAMATAKIS1=$(awk '{print $1}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy.stamatakis)

for j in {1..10}
do
sbatch -c $NC --output=$PHYL_DIR/logFiles/raxml-phylo_maxmiss$i.pars$j.oe $SCRIPTS_DIR/raxml-phylo_pars.sh $NC $OUTGROUP $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.pars$j.noinv.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4
sbatch -c $NC --output=$PHYL_DIR/logFiles/raxml-phylo_maxmiss$i.rand$j.oe $SCRIPTS_DIR/raxml-phylo_rand.sh $NC $OUTGROUP $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $PHYL_DIR/maxmiss$i/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.rand$j.noinv.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4

done


#2 for others (0.5, 0.75, 0.95) with checkpointing
MAXMISS=0.75
NJOBS=5
START_TREE=rand

STAMATAKIS1=$(awk '{print $1}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$MAXMISS.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$MAXMISS.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$MAXMISS.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$MAXMISS.noinv.phy.stamatakis)

for j in {1..10}
do

OUT_DIR=$PHYL_DIR/maxmiss$MAXMISS/$START_TREE/run$j
mkdir -p $OUT_DIR

for i in $(seq 1 $NJOBS)
do
#logging
echo "submitting job $i"
#special case: if it's the first script
if [ $i == 1 ]
then
#declare dir for saving checkpoint files to (it will be created in gphocs script)
CHECK_OUT=$OUT_DIR/checkpoint_$i
#submit job and save printed text to variable
SUBTEXT=$(sbatch -c $NC --output=$PHYL_DIR/logFiles/raxml-phylo_maxmiss$MAXMISS.$START_TREE.run$j.job$i.oe $SCRIPTS_DIR/raxml_checkpoint_$START_TREE.sh $NC $OUTGROUP $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$MAXMISS.noinv.phy $PHYL_DIR/maxmiss$MAXMISS/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$MAXMISS.$START_TREE.run$j.noinv.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4 $CHECK_OUT)
#extract sbatch ID and store in variable
declare RUNID_${i}=${SUBTEXT##* }

#otherwise: submit rerun scripts
else
#assign input dir (which is output dir of previous iteration)
CHECK_IN=$CHECK_OUT
#declare dir for saving checkpoint files to (it will be created in gphocs script)
CHECK_OUT=$OUT_DIR/checkpoint_$i
#get name of variable that stored sbatch ID in previous iteration of the loop
VARNAME=RUNID_$(( ${i} - 1 ))
#submit rerun and save again printed text to variable
echo ${!VARNAME}
SUBTEXT=$(sbatch -c $NC --output=$PHYL_DIR/logFiles/raxml-phylo_maxmiss$MAXMISS.$START_TREE.run$j.job$i.oe --dependency=afterany:${!VARNAME} $SCRIPTS_DIR/raxml_checkpoint_rerun.sh $CHECK_IN $CHECK_OUT $NC)
#see above
declare RUNID_${i}=${SUBTEXT##* }
fi
done
done


#3) for each maxmiss0.25-0.95 chose best scoring likelihood tree and do 100 bootstrap replicates in parallel















#########################
#same pipeline but only for focal individuals:
IND_FILE=/scratch/projects/nib00015/phylogenomicsProject/vcfFiltering/final/reducedSet.txt
bcftools view -S $IND_FILE $VCF_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.vcf > $VCF_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.vcf
sbatch --account=nib00015 --output=$PHYL_DIR/logFiles/prepareForRaxml_filtered04-hard_threshold0.25_reducedSet.oe $SCRIPTS_DIR/prepareForRaxml_SLURM.sh $VCF_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.vcf $PHYL_DIR
#CONTINUE HERE and check correctness of scripts#
NC=60
OUTGROUP="Mirzazaz_DLC2316m,Mirzazaz_DLC319m"
mkdir -p $PHYL_DIR/raxml/filtered04-hard_threshold0.25_reducedSet
STAMATAKIS1=$(awk '{print $1}' $PHYL_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $PHYL_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $PHYL_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $PHYL_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.noinv.phy.stamatakis)
sbatch --account=nib00015 -c $NC --output=$PHYL_DIR/logFiles/$PREFIX.filtered04-hard.raxml-phylo_0.25_reducedSet.oe $SCRIPTS_DIR/raxml-phylo.sh $NC $OUTGROUP $PHYL_DIR/$PREFIX.filtered04-hard.filtered06.threshold0.25.reducedSet.noinv.phy $PHYL_DIR/raxml/filtered04-hard_threshold0.25_reducedSet/$PREFIX.filtered04-hard.threshold0.25.reducedSet.noinv.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4
##########################


##compiling raxml-ng (not necessary because one can download precompiled versions)
#first install updated version of cmake locally then run commands specified on raxml-ng installation instructions (cmake needs to be run with the extra -DCMAKE_CXX_COMPILER=g++ after loading newer gcc module (e.g., gcc8))


###################################################
######Running SVDQ for SNPs with inds as tips######
###################################################

for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7 #0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.3 0.35 0.4
do
sbatch --output=$PHYL_DIR/logFiles/prepareForSVDQ_maxmiss$i.oe $SCRIPTS_DIR/prepareForSVDQ_SLURM.sh $VCF_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.thinned.vcf $PHYL_DIR
done


##creating locus partition block file not necessary since we run SVDquartets on SNPs


##create taxon partitions block file (is the same for each threshold) (for assignments to populations, create file manually)
echo "BEGIN SETS;" > $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
echo -e "\t TAXPARTITION SPECIES =" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt

for i in $(awk '{ print $1 }' $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss0.05.phy | tail -n+2)
do
echo -e "\t\t${i}:${i}," >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
done
#remove last comma
sed -i '$ s/,$//' $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt

echo -e "\t\t;" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
echo "END;" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
echo "" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt


#create paup block file
for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7 #0.05 0.25 0.5 0.75 0.95
do
SEED=$RANDOM
NQ=20000000
echo "BEGIN PAUP;" > $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt
echo -e "\toutgroup SPECIES.Cheirogaleusmaj_106245 SPECIES.Cheirogaleusmed_106354 SPECIES.Cheirogaleuscro_106189;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt
echo -e "\tset root=outgroup outroot=monophyl;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt
echo -e "\tsvdq nthreads=80 nquartets=${NQ} evalQuartets=random taxpartition=SPECIES bootstrap=standard seed=${SEED};" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt
echo -e "\tsavetrees format=Newick file=$PREFIX.maxmiss$i.thinned.svdq.tre savebootp=nodelabels;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt
echo -e "\tquit;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt
echo "END;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt

#concatenate nexus alignment file and partitions files (for individual partitions and for locus partitions)
cat $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.thinned.nex $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.txt > $PHYL_DIR/svdq/$PREFIX.concatenated.paup.maxmiss$i.thinned.nex
done


###run svdquartets in paup*
for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7 #0.05 0.25 0.5 0.75 0.95
do
NC=80
sbatch -c $NC --account=nib00015 --output=$PHYL_DIR/logFiles/$PREFIX.svdq.maxmiss$i.oe $SCRIPTS_DIR/svdq-phylo.sh $PHYL_DIR/svdq/$PREFIX.concatenated.paup.maxmiss$i.thinned.nex $PHYL_DIR/svdq/$PREFIX.concatenated.paup.maxmiss$i.thinned.nex.log
done


#######################################################
######Running SVDQ for SNPs with species as tips######
#######################################################

#create taxon partitions block file manually 
SPECIES_PARTITIONS=$PHYL_DIR/svdq/$PREFIX.taxonPartitions.species.txt


#create paup block file
for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7
do

SEED=$RANDOM
NQ=20000000
echo "BEGIN PAUP;" > $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt
echo -e "\toutgroup SPECIES.Cheirogaleusmaj_106245 SPECIES.Cheirogaleusmed_106354 SPECIES.Cheirogaleuscro_106189;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt
echo -e "\tset root=outgroup outroot=monophyl;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt
echo -e "\tsvdq nthreads=80 nquartets=${NQ} evalQuartets=random taxpartition=SPECIES bootstrap=standard seed=${SEED};" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt
echo -e "\tsavetrees format=Newick file=$PREFIX.maxmiss$i.thinned.speciesAssignments.svdq.tre savebootp=nodelabels;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt
echo -e "\tquit;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt
echo "END;" >> $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt

#concatenate nexus alignment file and partitions files (for individual partitions and for locus partitions)
cat $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.thinned.nex $PHYL_DIR/svdq/$PREFIX.taxonPartitions.species.txt $PHYL_DIR/svdq/$PREFIX.paup.maxmiss$i.speciesAssignments.txt > $PHYL_DIR/svdq/$PREFIX.concatenated.paup.maxmiss$i.thinned.speciesAssignments.nex
done

#run svdquartets in paup*
for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7 #0.05 0.25 0.5 0.75 0.95
do
NC=80
sbatch -c $NC --output=$PHYL_DIR/logFiles/$PREFIX.svdq.maxmiss$i.speciesAssignments.oe $SCRIPTS_DIR/svdq-phylo.sh $PHYL_DIR/svdq/$PREFIX.concatenated.paup.maxmiss$i.thinned.speciesAssignments.nex $PHYL_DIR/svdq/$PREFIX.concatenated.paup.maxmiss$i.thinned.speciesAssignments.nex.log
done


#####################################
#####running SplitsTree for SNPS#####
#####################################

#in order to run splitstree, nexus file has to be present
#run splitstree with Jelmer's script:

mkdir -p $PHYL_DIR/splitstree

for i in 0.05 0.25 0.5 0.75 0.95 #0.1 0.15 0.2 0.6 0.65 0.7 #0.05 0.25 0.5 0.75 0.95
do
#python $AMAS convert -i $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.phy -f phylip -u nexus -d dna
#mv $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.phy-out.nex $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.nex

OUT_FILE=$PHYL_DIR/splitstree/$PREFIX.splitstree.maxmiss$i.out.nex
sbatch --output=$PHYL_DIR/logFiles/$PREFIX.splitstree.maxmiss$i.oe $SCRIPTS_DIR/splitstree_run.sh $PHYL_DIR/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.nex $OUT_FILE 
done

#the outfile can then be opened in the GUI version of SplitsTree and a NeighborNet network can be created


##################################
######general locus processing####
##################################
SCRIPTS_DIR=/home/nibtve93/scripts/phylogeneticInference
PREFIX=final_15
LOCUS_DIR=/scratch/projects/nib00015/phylogenomicsProject/gatk/locusExtraction/final_reducedSet/fasta/reducedSet_bylocus_final_15
ALIGNMENT_DIR=/scratch/projects/nib00015/phylogenomicsProject/gatk/alignments/$PREFIX
PHYL_DIR=/scratch/projects/nib00015/phylogenomicsProject/gatk/phylogeneticInference/$PREFIX

mkdir -p $PHYL_DIR/logFiles
mkdir -p $PHYL_DIR/raxml
mkdir -p $PHYL_DIR/svdq
mkdir -p $ALIGNMENT_DIR/muscle
mkdir -p $ALIGNMENT_DIR/statistics
mkdir -p $ALIGNMENT_DIR/commands

#align loci with muscle
NC=32
sbatch -c $NC --output=$ALIGNMENT_DIR/muscle/muscle.oe $SCRIPTS_DIR/muscle.sh $LOCUS_DIR $ALIGNMENT_DIR $NC

#replace headers to remove locus IDs (otherwise concatenation will not work)
sbatch -c $NC --account=nib00015 --output=$ALIGNMENT_DIR/muscle/replaceHeaders.oe $SCRIPTS_DIR/replaceHeaders.sh $ALIGNMENT_DIR $NC

#calculate statistics for single locus alignments
$SCRIPTS_DIR/alignmentStats.sh $AMAS $ALIGNMENT_DIR/muscle $ALIGNMENT_DIR/statistics/muscleStatistics.txt

#concatenate alignment using AMAS
python $AMAS concat -i $ALIGNMENT_DIR/muscle/*.muscle.fa -t $ALIGNMENT_DIR/$PREFIX.concatenated.nex -p $ALIGNMENT_DIR/$PREFIX.partitions.txt -f fasta -d dna -u nexus -c 4
#replace N by ?
head -n 6 $ALIGNMENT_DIR/$PREFIX.concatenated.nex #check whether 6 lines is really the length of the header
awk 'NR>=7 {gsub("N","?",$2)}1' $ALIGNMENT_DIR/$PREFIX.concatenated.nex > $ALIGNMENT_DIR/$PREFIX.concatenated.nex.tmp
mv $ALIGNMENT_DIR/$PREFIX.concatenated.nex.tmp $ALIGNMENT_DIR/$PREFIX.concatenated.nex

##################################
######Running RAxML for loci######
##################################

#convert alignment from nexus to phylip format
cd $ALIGNMENT_DIR/
python $AMAS convert -i $ALIGNMENT_DIR/$PREFIX.concatenated.nex -f nexus -u phylip -d dna
mv $ALIGNMENT_DIR/$PREFIX.concatenated.nex-out.phy $ALIGNMENT_DIR/$PREFIX.concatenated.phy

#remove invariant sites (run over cluster because of high memory consumption)
sbatch -c 1 --output=$ALIGNMENT_DIR/removeInvariants.oe $SCRIPTS_DIR/removeInvariants.sh $ASCBIAS $ALIGNMENT_DIR/$PREFIX.concatenated.phy $ALIGNMENT_DIR/$PREFIX.concatenated.noinv.phy

#run raxml-ng
NC=80
OUTGROUP="Mirzazaz_DLC2316m,Mirzazaz_DLC319m"
STAMATAKIS1=$(awk '{print $1}' $ALIGNMENT_DIR/$PREFIX.concatenated.noinv.phy.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $ALIGNMENT_DIR/$PREFIX.concatenated.noinv.phy.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $ALIGNMENT_DIR/$PREFIX.concatenated.noinv.phy.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $ALIGNMENT_DIR/$PREFIX.concatenated.noinv.phy.stamatakis)
sbatch -c $NC --output=$PHYL_DIR/raxml/raxml-ng.oe $SCRIPTS_DIR/raxml-ng.sh $NC $OUTGROUP $ALIGNMENT_DIR/$PREFIX.concatenated.noinv.phy $PHYL_DIR/raxml/$PREFIX.concatenated.noinv.phy.tre $STAMATAKIS1 $STAMATAKIS2 $STAMATAKIS3 $STAMATAKIS4



##alternative estimations/models
#run raxml-ng without ascertainment bias correction and without partitioning (GTR+G model)
sbatch --account=nib00015 -c $NC --output=$PHYL_DIR/raxml/raxml-ng_vanilla.oe $SCRIPTS_DIR/raxml-ng_vanilla.sh $NC $OUTGROUP $ALIGNMENT_DIR/$PREFIX.concatenated.phy $PHYL_DIR/raxml/$PREFIX.concatenated.phy.vanilla.tre 100

#run raxml-ng without ascertainment bias correction but with partitioning (GTR+G model)
#takes a lot of memory, therefore runs on large40 instead of medium40
sed -e 's/^/GTR+G, /' $ALIGNMENT_DIR/$PREFIX.partitions.txt > $ALIGNMENT_DIR/$PREFIX.partitions.models.txt
sbatch --account=nib00015 -c $NC --output=$PHYL_DIR/raxml/raxml-ng_partitioned.oe $SCRIPTS_DIR/raxml-ng_partitioned.sh $NC $OUTGROUP $ALIGNMENT_DIR/$PREFIX.concatenated.phy $PHYL_DIR/raxml/$PREFIX.concatenated.phy.partitioned.tre 100 $ALIGNMENT_DIR/$PREFIX.partitions.models.txt

#run iq-tree with modelfinder without partitioning
sbatch --account=nib00015 -c $NC --output=$PHYL_DIR/iqtree/iqtree_vanilla.oe $SCRIPTS_DIR/iqtree_vanilla.sh $ALIGNMENT_DIR/$PREFIX.concatenated.phy $NC $PHYL_DIR/iqtree/vanilla

#run iq-tree with modelfinder and partitioning
mkdir -p $PHYL_DIR/iqtree
sed -e 's/^/DNA, /' $ALIGNMENT_DIR/$PREFIX.partitions.txt > $ALIGNMENT_DIR/$PREFIX.partitions.iqtree.txt
sbatch --account=nib00015 -c $NC --output=$PHYL_DIR/iqtree/iqtree_partitioned.oe $SCRIPTS_DIR/iqtree_partitioned.sh $ALIGNMENT_DIR/$PREFIX.concatenated.phy $NC $ALIGNMENT_DIR/$PREFIX.partitions.iqtree.txt $PHYL_DIR/iqtree/partitioned


#################################
######Running SVDQ for loci######
#################################

##create locus partitions block file
echo "BEGIN SETS;" > $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt
echo -e "\t CHARPARTITION LOCI =" >> $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt

while read line
do
locus_name=$(cut -f1 -d' ' <<< $line)
locus_coordinates=$(cut -f3 -d' ' <<< $line)
echo "Processing locus $locus_name"
echo -e "\t\t${locus_name}:${locus_coordinates}," >> $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt
done < $ALIGNMENT_DIR/$PREFIX.partitions.txt
#remove last comma
sed -i '$ s/,$//' $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt

echo -e "\t\t;" >> $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt
echo "END;" >> $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt
echo "" >> $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt

##create taxon partitions block file (for assignments to populations, create file manually)
echo "BEGIN SETS;" > $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
echo -e "\t TAXPARTITION SPECIES =" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt

for i in $(awk '{ print $1 }' $ALIGNMENT_DIR/$PREFIX.concatenated.phy | tail -n+2)
do
echo -e "\t\t${i}:${i}," >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
done
#remove last comma
sed -i '$ s/,$//' $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt

echo -e "\t\t;" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
echo "END;" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt
echo "" >> $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt

#create paup block file (adjust outgroups!)
#nthreads medium40 partition has 40 CPU cores on 2 sockets (which can be hyperthreaded to 80 logical cores)
SEED=$RANDOM
echo "BEGIN PAUP;" > $PHYL_DIR/svdq/$PREFIX.paup.txt
echo -e "\toutgroup SPECIES.Mirzazaz_DLC2316m SPECIES.Mirzazaz_DLC319m;" >> $PHYL_DIR/svdq/$PREFIX.paup.txt
echo -e "\tset root=outgroup outroot=monophyl;" >> $PHYL_DIR/svdq/$PREFIX.paup.txt
echo -e "\tsvdq nthreads=80 evalQuartets=all taxpartition=SPECIES loci=LOCI bootstrap=standard seed=${SEED};" >> $PHYL_DIR/svdq/$PREFIX.paup.txt
echo -e "\tsavetrees format=Newick file=$PREFIX.svdq.bootstrap.tre savebootp=nodelabels;" >> $PHYL_DIR/svdq/$PREFIX.paup.txt
echo -e "\tquit;" >> $PHYL_DIR/svdq/$PREFIX.paup.txt
echo "END;" >> $PHYL_DIR/svdq/$PREFIX.paup.txt

#concatenate nexus alignment file and partitions files (for individual partitions and for locus partitions)
cat $ALIGNMENT_DIR/$PREFIX.concatenated.nex $PHYL_DIR/svdq/$PREFIX.locusPartitions.txt $PHYL_DIR/svdq/$PREFIX.taxonPartitions.txt $PHYL_DIR/svdq/$PREFIX.paup.txt > $PHYL_DIR/svdq/$PREFIX.concatenated.paup.nex

###run svdquartets in paup*
NC=80
sbatch --account=nib00015 -c $NC --output=$PHYL_DIR/svdq/$PREFIX.svdq.oe $SCRIPTS_DIR/svdq.sh $PHYL_DIR/svdq/$PREFIX.concatenated.paup.nex $PHYL_DIR/svdq/$PREFIX.concatenated.paup.nex.log


