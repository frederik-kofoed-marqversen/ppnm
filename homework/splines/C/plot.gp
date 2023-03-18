set terminal svg

set output "Plot.svg"
set key top right
set title "Spline interpolation of Lorenz distribution from 10 samples"
set xlabel "x"
set ylabel "f(x)"
set tics out
plot "data.txt" index 0 with points title "Sample points" ,\
    "data.txt" index 1 using ($1):($3) with lines title "Cubic spline" ,\
    "data.txt" index 1 using ($1):($4) with lines title "Derivative" ,\
    "data.txt" index 1 using ($1):($5) with lines title "Integral"