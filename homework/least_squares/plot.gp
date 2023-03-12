set terminal svg

set output "plot.svg"
set key top right
set xlabel "time (days)"
set ylabel "Activity (relative units)"
set tics out
set title "Thorium X decay"
plot "Out.txt" index 0 with errorbars title "Experiment",\
     "Out.txt" index 1 using 1:2 with lines title "Lst.Sq. fit",\
     "Out.txt" index 2 using 1:2 with lines dashtype 2 title "Lst.Sq. bounds"