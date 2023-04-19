set terminal svg background "white"

set key top right
set xlabel "Number of samples (N)"
set ylabel "Error"
set tics out
set title "Error scaling"

set fit quiet
set fit logfile '/dev/null'
f(x) = a / sqrt(x)
fit f(x) 'scaling.txt' index 0 using 1:3 via a
g(x) = b / sqrt(x)
fit g(x) 'scaling.txt' index 0 using 1:(abs($2-13.0/20.0)) via b
h(x) = c / sqrt(x)
fit h(x) 'scaling.txt' index 1 using 1:(abs($2-13.0/20.0)) via c

plot [1e4:1e5] "scaling.txt" index 0 using 1:3 with points title "Plain Monte Carlo estimate", \
    "scaling.txt" index 0 using 1:(abs($2-13.0/20.0)) with points title "Plain Monte Carlo true", \
    f(x) with lines linecolor "black" notitle, \
    g(x) with lines linecolor "black" notitle, \
    "scaling.txt" index 1 using 1:(abs($2-13.0/20.0)) with points title "Low descrepancy true", \
    h(x) with lines linecolor "black" title "Fit to: f(N)=a/sqrt(N)"
