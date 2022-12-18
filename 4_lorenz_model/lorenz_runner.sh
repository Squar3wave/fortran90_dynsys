#!/bin/bash

i=0

if [ ! -d "data"  ];
then
  mkdir data
fi

file="parameters.txt"

# Lorenz parameters

s=$(sed -n '2p' $file)
r=$(sed -n '3p' $file)
b=$(sed -n '4p' $file)

# pos_1 (in order x, y and z) for R-K 4

ps_1_0=$(sed -n '7p' $file)
ps_1_1=$(sed -n '8p' $file)
ps_1_2=$(sed -n '9p' $file)

# pos_2 (in order x, y and z) for single Lyapunov coeff

ps_2_0=$(sed -n '12p' $file)
ps_2_1=$(sed -n '13p' $file)
ps_2_2=$(sed -n '14p' $file)

# dpos_2 (in order x, y and z) for single Lyapunov coeff

dp_2_0=$(sed -n '17p' $file)
dp_2_1=$(sed -n '18p' $file)
dp_2_2=$(sed -n '19p' $file)

# pos_a (in order x, y and z) for triple Lyapunov coeff (same for a,b and c)

ps_a_0=$(sed -n '22p' $file)
ps_a_1=$(sed -n '23p' $file)
ps_a_2=$(sed -n '24p' $file)

# dpos_a (in order x, y and z) for triple Lyapunov coeff

dp_a_0=$(sed -n '27p' $file)
dp_a_1=$(sed -n '28p' $file)
dp_a_2=$(sed -n '29p' $file)

# dpos_b (in order x, y and z) for triple Lyapunov coeff

dp_b_0=$(sed -n '32p' $file)
dp_b_1=$(sed -n '33p' $file)
dp_b_2=$(sed -n '34p' $file)

# dpos_c (in order x, y and z) for triple Lyapunov coeff

dp_c_0=$(sed -n '37p' $file)
dp_c_1=$(sed -n '38p' $file)
dp_c_2=$(sed -n '39p' $file)

d=42

for i in {0..12}
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

# program     s  r  b  x       y       z

./mod_lorenz $s $r $b $ps_1_0 $ps_1_1 $ps_1_2\
                      $ps_2_0 $ps_2_1 $ps_2_2\
                      $dp_2_0 $dp_2_1 $dp_2_2\
                      $ps_a_0 $ps_a_1 $ps_a_2\
                      $dp_a_0 $dp_a_1 $dp_a_2\
                      $ps_a_0 $ps_a_1 $ps_a_2\
                      $dp_b_0 $dp_b_1 $dp_b_2\
                      $ps_a_0 $ps_a_1 $ps_a_2\
                      $dp_c_0 $dp_c_1 $dp_c_2\
             $dist_0 $dist_1 $dist_2 $dist_3 $dist_4\
             $dist_5 $dist_6 $dist_7 $dist_8 $dist_9 $dist_10 $dist_11 $dist_12
             
echo ""
echo "========================================================================="
echo ""
echo "Gnuplot"
echo ""
echo "========================================================================="
echo ""

gnuplot lorenz_plotter.plt
