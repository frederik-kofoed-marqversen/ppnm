set terminal svg dynamic size 600,900 

set key bottom right
set ylabel "Energy error (Hartree)"
set format y "%.1e"
set multiplot layout 3,1
    set xlabel "r_{min} (Bohr)"
    plot [0.10:0] "convergence.txt" index 0 using 3:($1+0.5) with lines notitle
    set xlabel "r_{max} (Bohr)"
    plot "convergence.txt" index 1 using 2:($1+0.5) with lines notitle
    set xlabel "absolute precision of ODE-solver"
    set log x 10
    set format x "%.0e"
    plot [1e-2:1e-4] "convergence.txt" index 2 using 4:($1+0.5) with lines notitle