set terminal svg

set output "Out/states.svg"
set key top right
set xlabel "r (Bohr)"
set ylabel "f(r)"
set tics out
set title "3 lowest Hydrogen eigenstates"
plot for [idx=0:2] "states.data" index idx with lines title "state: ".idx

set output "Out/dr_conv.svg"
set key top left
set xlabel "dr (Bohr)"
set ylabel "energy"
set tics out
set title "Convergence in stepsize"
plot [0.45: 0][-0.51: -0.47] "dr.data" using ($3):($1) with points

set output "Out/r_max_conv.svg"
set key top right
set xlabel "r_{max} (Bohr)"
set ylabel "energy"
set tics out
set title "Convergence in r_{max}"
plot [1: 10][-0.6: 0] "r_max.data"  using ($2):($1) with points