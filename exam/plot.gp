set terminal svg background "white"

set key top left
set fit quiet
set fit logfile '/dev/null'
f(x) = a * x**3
g(x) = b * x**3
fit f(x) 'times.txt' index 0 via a
fit g(x) 'times.txt' index 1 via b
plot [0:1000] "times.txt" index 0 with points linecolor "black" title "Hessenberg decomposition",\
    "times.txt" index 1 with points linecolor "red" title "QR factorisation",\
    f(x) with lines linecolor rgb "black" title "fit to: a*x^3",\
    g(x) with lines linecolor "red" title "fit to: a*x^3"
