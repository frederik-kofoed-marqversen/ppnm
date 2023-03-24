set terminal svg

set fit quiet
set fit logfile '/dev/null'
f(x) = a * x**3
fit f(x) 'times.data' via a
plot "times.data" with points,\
    f(x) with lines title "fit to: a*x^3"