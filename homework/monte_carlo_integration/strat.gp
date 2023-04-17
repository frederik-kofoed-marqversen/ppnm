set terminal svg background "white"

set title "5000 stratified samples"
set tics out
set size ratio -1

plot "strat.txt" with points pointtype 6 linecolor "red" notitle, \
    sqrt(1.0**2 - x**2) with lines linewidth 2 dashtype 2 linecolor "black" notitle, \
    sqrt(0.5**2 - x**2) with lines linewidth 2 dashtype 2 linecolor "black" notitle