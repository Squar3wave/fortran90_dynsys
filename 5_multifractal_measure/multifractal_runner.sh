#!/bin/bash

if [ ! -d "data" ];
then
  mkdir data
fi

echo "========================================================================="
echo ""
echo "Fortran program"
echo ""
echo "========================================================================="
echo ""

./multifract_albero


echo ""
echo "========================================================================="
echo ""
echo "Gnuplot"
echo ""
echo "========================================================================="
echo ""

gnuplot multifractal_plotter.plt

