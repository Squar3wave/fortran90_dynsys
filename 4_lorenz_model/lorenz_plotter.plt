set term png
set out "lorenz.png"


sp 'data/rk4_lorenz.dat' w l, 'data/stz.dat'


unset term
