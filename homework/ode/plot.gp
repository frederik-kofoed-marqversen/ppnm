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
set samples 101
plot "data.txt" index 0 using 1:2 smooth csplines with lines linewidth 2 title "theta(t)", \
     "data.txt" index 0 using 1:3 smooth csplines with lines linewidth 2 title "omega(t)"

set output "Lotka_Volterra.svg"
set key top left
set xlabel "t"
set tics out
set title "Lotka-Volterra system"
unset xzeroaxis
unset grid
set samples 500
plot "data.txt" index 1 using 1:2 smooth csplines with lines linewidth 2 title "x(t)", \
     "data.txt" index 1 using 1:3 smooth csplines with lines linewidth 2 title "y(t)"

set output "Three_body.svg"
set key top left
unset xlabel
unset ylabel
set tics out
set title "Three body system"
set xzeroaxis
set yzeroaxis
set grid
set size ratio -1
plot [-1.5:1.5][-1:1] "data.txt" index 2 using 1:2 with lines linewidth 2 title "r_1(t)", \
     "data.txt" index 2 using 3:4 with lines linewidth 2 title "r_2(t)", \
     "data.txt" index 2 using 5:6 with lines linewidth 2 title "r_3(t)"