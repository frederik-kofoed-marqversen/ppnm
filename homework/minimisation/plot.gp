set terminal svg

set key top right
set xlabel "E[GeV]"
set ylabel "signal σ(E)"
set tics out
set title "Hydrogen ground state"
plot "higgs.data" using 1:2 every ::1 with points title "Higgs data", \
        "Out.txt" index 1 with lines title "Fit to Breit-Wigner"