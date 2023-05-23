set terminal svg background "white"

set output "Plot.svg"
set key top right
set title "Spline interpolation from 20 samples"
set xlabel "x"
set ylabel "f(x)"
set tics out
plot "data.txt" using ($1):($2) with lines title "spline: sin(x)",\
     "data.txt" using ($1):($4) with lines title "integral: 1-cos(x)",\
     "data.txt" using ($1):($3) with lines title "derivative: cos(x)"