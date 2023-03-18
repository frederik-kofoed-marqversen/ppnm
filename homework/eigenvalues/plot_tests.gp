set terminal svg

size = system("tail -n 1 times.txt | cut -f 1")
time = system("tail -n 1 times.txt | cut -f 2")
plot "times.txt" with points,\
    time/(size**3)*x**3 with lines title "f(x) = a x^3"