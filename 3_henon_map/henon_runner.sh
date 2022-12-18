#!/bin/bash

if [ ! -d "data" ];
then
  mkdir data
fi

file="parameters.txt"

a=$(sed -n '2p' $file)
b=$(sed -n '3p' $file)

x0_0=$(sed -n '6p' $file)
x0_1=$(sed -n '7p' $file)

x1_0=$(sed -n '10p' $file)
x1_1=$(sed -n '11p' $file)

dx1_0=$(sed -n '14p' $file)
dx1_1=$(sed -n '15p' $file)

x2_0=$(sed -n '18p' $file)
x2_1=$(sed -n '19p' $file)

dx2_0=$(sed -n '22p' $file)
dx2_1=$(sed -n '23p' $file)

dx3_0=$(sed -n '26p' $file)
dx3_1=$(sed -n '27p' $file)

d=30
 
 for i in {0..10}
 do
     l=$(($d+$i))
     declare dist_$i=$(sed -n "$l p" $file)
 done


echo "========================================================================="
echo ""
echo "Fortran program"
echo ""
echo "========================================================================="
echo ""

./map_henon $a $b $x0_0 $x0_1\
                  $x1_0 $x1_1 $dx1_0 $dx1_1\
                  $x2_0 $x2_1 $dx2_0 $dx2_1\
                  $x2_0 $x2_1 $dx3_0 $dx3_1\
            $dist_0 $dist_1 $dist_2 $dist_3 $dist_4\
            $dist_5 $dist_6 $dist_7 $dist_8 $dist_9 $dist_10 $dist_11 

echo ""
echo "========================================================================="
echo ""
echo "Gnuplot"
echo ""
echo "========================================================================="
echo ""

gnuplot henon_plotter.plt
