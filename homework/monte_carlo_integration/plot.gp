set terminal svg

set key top right
set xlabel "Number of samples (N)"
set ylabel "Error"
set tics out
set title "Error scaling"

set fit quiet
set fit logfile '/dev/null'
f(x) = a / sqrt(x)
fit f(x) 'scaling.data' index 0 via a

plot "scaling.data" index 0 with points title "Plain Monte Carlo", \
    f(x) with lines title "Fit to: f(N)=a/sqrt(N)"
