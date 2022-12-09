#!/bin/bash

if [ ! -d "data" ];
then
  mkdir data
fi

./lyap

gnuplot logistic_plotter.plt
