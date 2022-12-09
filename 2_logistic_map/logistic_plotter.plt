set term png
set out "lyapunov_all.png"

set grid
set ylabel "lyapunov coefficients"
set xlabel 'time'
FILES = system('ls -1 data/data*.dat')
plot for [data in FILES] data with lines

unset term

set term png
set out "lyapunov_behaviour.png"

set grid
set ylabel "lyapunov coefficient"
set xlabel 'mu'
plot 'data/lyap_mu.dat' with lines

unset term
