set terminal svg background "white"

set key top right
set xlabel "r (Bohr)"
set ylabel "f(r)"
set tics out
set title "Hydrogen ground state"
plot    [0:14] "/dev/stdin" index 1 with points title "Numeric solution", \
        "/dev/stdin" index 1 smooth csplines with lines title "Ground state spline", \
        "/dev/stdin" index 2 smooth csplines with lines title "Exited state Spline", \
        "/dev/stdin" index 3 smooth csplines with lines title "Exited state Spline", \
        x*exp(-x) with lines linewidth 2 dashtype 2 linecolor "black" title "Analytical result"