#!/bin/bash


in=$1

total=`wc -l $in | awk '{ print $1 }'`

echo "Total elements in the database: " $total

active=`grep ligand $in | wc -l | awk '{ print $1 }'`

echo "Active elements in the database: " $active

selected1=`awk -v a=$total 'BEGIN { res = sprintf("%.0f", a/100); print res }'`
selected5=`awk -v a=$total 'BEGIN { res = sprintf("%.0f", a/20); print res }'`
selected10=`awk -v a=$total 'BEGIN { res = sprintf("%.0f", a/10); print res }'`


echo "Selected elements in 1%, 5% and 10% database: " $selected1 " " $selected5 " " $selected10


TP1=`head -$selected1 $in | grep ligand | wc -l | awk '{ print $1 }'`
TP5=`head -$selected5 $in | grep ligand | wc -l | awk '{ print $1 }'`
TP10=`head -$selected10 $in | grep ligand | wc -l | awk '{ print $1 }'`

echo "True positives in 1%, 5% and 10% database: " $TP1 $TP5 $TP10

EF1=`awk -v tp1=$TP1 -v sel1=$selected1 -v act=$active -v tot=$total 'BEGIN { res = sprintf("%.2f", (tp1/sel1)/(act/tot)); print res }'`
EF5=`awk -v tp5=$TP5 -v sel5=$selected5 -v act=$active -v tot=$total 'BEGIN { res = sprintf("%.2f", (tp5/sel5)/(act/tot)); print res }'`
EF10=`awk -v tp10=$TP10 -v sel10=$selected10 -v act=$active -v tot=$total 'BEGIN { res = sprintf("%.2f", (tp10/sel10)/(act/tot)); print res }'`


echo $EF1 " " $EF5 " " $EF10
