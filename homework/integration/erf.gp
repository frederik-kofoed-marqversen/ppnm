set terminal svg

set key top left
set xlabel "x"
set ylabel "erf(x)"
set tics out
set title "Error function"
set xzeroaxis
set yzeroaxis
set grid

plot [-3: 3][-1.1: 1.1] "erf.data" with lines