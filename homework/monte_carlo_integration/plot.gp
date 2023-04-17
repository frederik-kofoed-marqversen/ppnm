set terminal svg

set key top right
set xlabel "Number of samples (N)"
set ylabel "Error"
set tics out
set title "Error scaling"

set fit quiet
set fit logfile '/dev/null'
f(x) = a / sqrt(x)
fit f(x) 'scaling.txt' index 0 via a
g(x) = b / sqrt(x)
fit g(x) 'scaling.txt' index 1 via b

plot "scaling.txt" index 0 with points title "Plain Monte Carlo", \
    f(x) with lines title "Fit to: f(N)=a/sqrt(N)", \
    "scaling.txt" index 1 with points title "Low descrepancy", \
    g(x) with lines title "Fit to: f(N)=a/sqrt(N)"
