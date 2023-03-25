set terminal svg

set key top left
set fit quiet
set fit logfile '/dev/null'
f(x) = a * x**3
g(x) = b * x**3
fit f(x) 'times.data' index 0 via a
fit g(x) 'times.data' index 1 via b
plot "times.data" index 0 with points title "Jacobi cyclic",\
    "times.data" index 1 with points title "Optimised alg",\
    f(x) with lines title "fit to: a*x^3",\
    g(x) with lines title "fit to: a*x^3"