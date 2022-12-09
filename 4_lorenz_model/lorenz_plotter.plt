set term png
set out "lorenz.png"


sp 'rk4_lorenz.dat' w l, 'stz.dat'


unset term
