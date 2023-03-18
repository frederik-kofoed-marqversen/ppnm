set terminal svg background "white"
set key top left

size = system("tail -n 1 Out/times.data | cut -f 1")
time = system("tail -n 1 Out/times.data | cut -f 2")
plot [0:2100] "Out/times.data" with points \
            , time/(size**3)*x**3 with lines title "f(x) = a x^3"
