set terminal svg background "white"
set key top left

plot [0:2100] "out/times.data" with points \
            , 10.24/(2000**3)*x**3 with lines
