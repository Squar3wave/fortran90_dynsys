set term png
set out "final_system_ascending_sigma.png"

set grid 
set xlabel 'time'
set ylabel 'position'
p 'data_final/data_final_p.dat' w l

unset term


set term png
set out "final_system_descending_sigma.png"

set grid 
set xlabel 'time'
set ylabel 'position'
p 'data_final/data_final_m.dat' w l

unset term


set term png
set out "maxes_ascending_sigma.png"

set grid 
set xlabel 'sigma'
set ylabel 'maxes'
p 'maxes/max_p.dat' w l

unset term


set term png
set out "maxes_descending_sigma.png"

set grid 
set xlabel 'sigma'
set ylabel 'maxes'
p 'maxes/max_m.dat' w l

unset term


set term png
set out "maxes_compared.png"

set grid 
set xlabel 'sigma'
set ylabel 'maxes'
p 'maxes/max_m.dat' w l, 'maxes/max_p.dat' w l

unset term

