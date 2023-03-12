#!/bin/bash
cat > plot.gp << %EOF%
set terminal svg background "white"
set key top left

plot [0:2100] "out/times.data" with points \\
            , $1/(2000**3)*x**3 with lines
%EOF%