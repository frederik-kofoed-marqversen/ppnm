set terminal svg dynamic size 1200,900 background "white"

set key top right
set ylabel "|E_{calc} - E_{true}|   [Hartree]"
set format y "%.1e"
set logscale y 10
set multiplot layout 3,2 columnsfirst title "Energy convergence part B and C" font ",24"
# PART B
     set xlabel "r_{max}   [Bohr]"
     plot "/dev/stdin" index 1 using 2:(abs($1+0.5)) with lines title "Ground state - part B"
     
     set xlabel "r_{min}   [Bohr]"
     plot [0.08:0] "/dev/stdin" index 0 using 3:(abs($1+0.5)) with lines title "Ground state - part B"
    
     set xlabel "absolute precision of ODE-solver"
     set logscale x 10
     set format x "%.0e"
     plot [10**(-2.5):1e-4] "/dev/stdin" index 2 using 4:(abs($1+0.5)) with lines title "Ground state - part B"
     unset logscale
     unset format
# PART C
     set xlabel "r_{max}   [Bohr]"
     set format y "%.1e"
     set logscale y 10
     plot [5:30] "/dev/stdin" index 3 using 2:(abs($1+0.5)) with lines title "Ground state - part C"
     plot [5:30] "/dev/stdin" index 4 using 2:(abs($1+0.5/4)) with lines title "1^{st} excited state - part C"
     plot [5:30] "/dev/stdin" index 5 using 2:(abs($1+0.5/9)) with lines title "2^{nd} excited state - part C"