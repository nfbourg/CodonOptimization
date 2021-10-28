#!/bin/bash

input=$1
prefix=$2

ws=200
step=10
dir=/grid/home/nbourgeois/codon_optimization/DeepPASTA_gpu
dir_ss=$dir/generating_secondary_structure_from_sequence
dir_pas=$dir/polyA_site_prediction

len=`infoseq $input | awk '{print $6}' | tail -n +2`
name=`grep ">" $input | sed 's/>//'`

#>>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/software/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/software/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/data/software/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data/software/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
#<<< conda initialize <<<

conda activate DeepPASTA_test

# Preparing files
rm -f $prefix.fai

# insert1
rm -f $prefix.bed
rm -f $prefix.stag.bed

while read -r line;
do
    if [[ $line == \>* ]] ;
    then
        echo -n "$line" | sed 's/>//'  >> $prefix.bed
        echo -n "$line" | sed 's/>//'  >> $prefix.stag.bed
    else
        line=`echo $line | sed -e 's/^[[:space:]]*//'`
        len=${#line}
        strt=$(($len - $ws))
        echo " 1 $len" | sed 's/ /\t/g'  >> $prefix.bed
        echo " $strt $len" | sed 's/ /\t/g'  >> $prefix.stag.bed
    fi
done < "$input"
#end insert1

bedtools makewindows -b $prefix.bed -w $ws -s $step | awk '($3 - $2 + 1 ) == (ws+1) {print $0}' ws=$ws >> $prefix.stag.bed

bedtools getfasta -fi $input -bed $prefix.stag.bed -fo $prefix.stag.fasta

# Calculating secondary structure
input_stag=$prefix.stag.fasta

$dir_ss/RNAshapes \
-f  $input_stag -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > $prefix.ss.txt

python $dir_ss/combining_substructure.py -i $prefix.ss.txt -o $prefix.ss.comb.txt
python $dir_ss/filtering_number_of_ss.py -n 3 -i $prefix.ss.comb.txt -o $prefix.ss.filt.txt
python $dir_ss/shape_assign_per_nucleotide.py -c 3 -i $prefix.ss.filt.txt -o $prefix.ss.per_nt.txt

# Predicting Polyadenilation site (PAS)
input_ss=$prefix.ss.per_nt.txt

hdf5_link=DeepPASTA_polyA_site_learned.hdf5
hdf5_dest=$dir_pas/DeepPASTA_polyA_site_learned.hdf5


if [ -L ${hdf5_link} ] ; then
   if [ -e ${hdf5_link} ] ; then
	# Good link, do nothing
	echo "We already have a link to hd5"
   else
	# Broken link
	rm -f $hdf5_link
	ln -s $hdf5_dest
   fi
elif [ -e ${my_link} ] ; then
   # Not a link, do nothing
   echo "It seems that we have copy of hd5 in the current directory"
else
   # Missing
   rm -f $hdf5_link
   ln -s $hdf5_dest
fi

python $dir/DeepPASTA_polyA_site_prediction_testing.py \
-testSeq $input_stag \
-testSS $input_ss \
-o $prefix.report.txt

# Saving result
hiscore=`sort -k2,2gr $prefix.report.txt | head -n 1`
echo $hiscore > $prefix.pas.txt

