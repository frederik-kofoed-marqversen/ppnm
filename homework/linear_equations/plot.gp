set terminal svg background "white"
set key top left

set fit quiet
set fit logfile '/dev/null'
f(x) = a * x**3
fit f(x) 'Out/times.data' via a
plot [0:2100] "Out/times.data" with points, \
            f(x) with lines title "fit to: a*x^3"
