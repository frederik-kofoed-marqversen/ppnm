set terminal svg

set output "Plot.svg"
set key top right
set xlabel "time (days)"
set ylabel "Activity (relative units)"
set tics out
set title "Thorium X decay"
plot "out.txt" index 0 with errorbars title "Experiment",\
     "out.txt" index 1 using 1:2 with lines title "Lst.Sq. fit",\
     "out.txt" index 2 using 1:2 with lines dashtype 2 title "Lst.Sq. bounds"