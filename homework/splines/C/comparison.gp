set terminal svg

set output "Gnuplot_comparison.svg"
set key top right
set title "Spline interpolation from 10 samples"
set xlabel "x"
set ylabel "f(x)"
set tics out
plot "data.txt" index 1 using ($1):($3) with lines linecolor "black" title "My cubic spline" ,\
    "data.txt" index 0 smooth csplines dashtype 2 linewidth 3 linecolor "red" title "Gnuplot cubic spline" ,\
    "data.txt" index 0 with points pointtype 5 pointsize 0.8 linecolor "black" title "Sample points" ,\