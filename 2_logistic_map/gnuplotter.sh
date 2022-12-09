#!/bin/bash

gnuplot -p -e  "set grid; FILES = system('ls -1 data*.dat'); plot for [data in FILES] data with lines"
gnuplot -p -e  "set grid; set xlabel 'mu'; set ylabel 'lyapunov coefficient' ; plot 'lyap_mu.dat' with lines"
