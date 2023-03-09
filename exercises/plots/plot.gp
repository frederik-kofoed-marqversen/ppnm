set terminal svg

set output "out/gamma.svg"
set key bottom right
set xlabel "x"
set ylabel "y"
set tics out
set title "Factorial 'function'"
set xzeroaxis
set yzeroaxis
plot [-5: 5][-5: 5] \
    "out/factorials.data" with points pointtype 2 title "factorial", \
    "out/gamma.data" using ($1-1):($2) with lines linetype 1 title "shifted gamma function"

set output "out/lngamma.svg"
set key top left
set xlabel "x"
set ylabel "y"
set tics out
set title "Logarithmic Gamma function"
set xzeroaxis
set yzeroaxis
plot "out/lngamma.data" with lines linetype 1 title "y = ln(gamma(x))"