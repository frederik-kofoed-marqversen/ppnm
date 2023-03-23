set terminal svg

set output "Pendulum.svg"
set key bottom right
set xlabel "t"
# set ylabel "y"
set tics out
set title "Damped pendulum"
set xzeroaxis
set yzeroaxis
set grid
plot "data.txt" index 0 using ($1):($2) with lines linewidth 2 title "theta(t)", \
     "data.txt" index 0 using ($1):($3) with lines linewidth 2 title "omega(t)"

set output "Lotka_Volterra.svg"
set key top left
set xlabel "t"
set tics out
set title "Lotka-Volterra System"
unset xzeroaxis
unset grid
plot "data.txt" index 1 using ($1):($2) with lines linewidth 2 title "x(t)", \
     "data.txt" index 1 using ($1):($3) with lines linewidth 2 title "y(t)"