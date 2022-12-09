set term png
set out "multifractal_signal.png"

set title 'Multifractal signal'
p 'signal.dat' w l

unset out

set term png
set out "multifractal_moments.png"

set title 'Linear fitting for momentums in logarithmic scale'
set yrange [-1:]
f(x) = m1*x + q1
g(x) = m2*x + q2
h(x) = m3*x + q3
fit f(x) 'mft_lnm_1.dat' using 1:2 via m1,q1
fit g(x) 'mft_lnm_2.dat' using 1:2 via m2,q2
fit h(x) 'mft_lnm_3.dat' using 1:2 via m3,q3
p 'mft_lnm_1.dat', f(x), 'mft_lnm_2.dat', g(x), 'mft_lnm_3.dat', h(x)

unset out
