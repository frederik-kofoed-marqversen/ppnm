set terminal svg

set output "target/gamma.svg"
set key bottom right
set xlabel "x"
set ylabel "y"
set tics out
set title "Factorial 'function'"
set xzeroaxis
set yzeroaxis
plot [-5: 5][-5: 5] \
    "target/factorials.data" with points pointtype 2 title "factorial", \
    "target/gamma.data" using ($1-1):($2) with lines linetype 1 title "shifted gamma function"

set output "target/lngamma.svg"
set key top left
set xlabel "x"
set ylabel "y"
set tics out
set title "Logarithmic Gamma function"
set xzeroaxis
set yzeroaxis
plot "target/lngamma.data" with lines linetype 1 title "y = ln(gamma(x))"