set terminal svg

set key top right
set xlabel "r (Bohr)"
set ylabel "f(r)"
set tics out
set title "Hydrogen ground state"
plot "Out.txt" index 1 every ::1 with points title "Numeric solution", \
        "Out.txt" index 1 every ::1 smooth csplines with lines title "Spline", \
        x*exp(-x) with lines linewidth 2 dashtype 2 linecolor "black" title "Analytical result"