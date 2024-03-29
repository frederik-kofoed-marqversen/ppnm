set terminal svg background "white"

set key top right
set xlabel "E[GeV]"
set ylabel "signal σ(E)"
set tics out
set title "Hydrogen ground state"
plot "higgs.data" using 1:2:3 with yerrorbars title "Higgs data", \
        "Out.txt" index 1 with lines title "Fit to Breit-Wigner"