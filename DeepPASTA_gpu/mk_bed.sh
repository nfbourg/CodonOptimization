#!/bin/bash

input=$1
prefix=$2
ws=200

rm -f $prefix.bed

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
