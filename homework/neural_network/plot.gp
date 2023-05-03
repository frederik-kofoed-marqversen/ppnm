set terminal svg background "white" 

# Part A
set output "Curve_fit.svg"
set key top left
set tics out
set title "Neural network curve fit with 4 hidden neurons"
plot [-1:1] cos(5.0*x - 1.0) * exp(-x*x) with lines linewidth 2 dashtype 2 linecolor "black" title "Analytical function", \
    "data.txt" index 0 with points title "Training data", \
    "data.txt" index 1 with lines title "Neural network fit"

# Part C (and implicitly B)
set output "Differential_equation.svg"
set key top left
set tics out
set title "Neural network solution to differential equation with 6 hidden neurons"
plot [-2*pi:2*pi] cos(x) with lines linewidth 2 dashtype 2 linecolor "black" title "Analytical solution: cos(x)", \
    "data.txt" index 2 with lines title "Neural network"