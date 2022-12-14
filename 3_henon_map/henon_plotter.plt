set term png
set out "henon_euler.png"

set title "Henon map"
p 'data/euler_henon_1.dat'

unset term

 
set term png
set out "henon_linfit.png"

f(x) = m*x+q
fit f(x) 'data/gp_henon_log.dat' using 1:2 via m,q

set title "linear fit"
p f(x), 'data/gp_henon_log.dat' 

unset term
