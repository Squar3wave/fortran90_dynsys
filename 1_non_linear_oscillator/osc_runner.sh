#!/bin/bash

if [ ! -d "maxes" ];
then
  mkdir maxes
fi

if [ ! -d "data_final" ];
then
  mkdir data_final
fi

if [ ! -d "data_m" ];
then
  mkdir data_m
fi

if [ ! -d "data_p" ];
then
  mkdir data_p
fi

./osc_non_lin

gnuplot oscillator_plotter.plt
