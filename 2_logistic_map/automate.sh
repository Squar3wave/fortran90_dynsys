#!/bin/bash

echo "start cleaning"
rm data_*
echo "done cleaning, start compiling"
make
echo "done compiling, starting program"
./lyap
./gnuplotter.sh
echo "done!"
