#set terminal latex
#set terminal wxt
set terminal epslatex size 15cm,8cm color colortext
set output "testCase2-cputvscellc-s-x-T-2000-t-10-viscrat-".muw."_".muo.".tex"
set format xy "$%g$"
#set format xy "%g"
#set title "CPU time vs time step length"
set xlabel "Cell count"
set ylabel "CPU time [seconds]"
set ylabel rotate left
set logscale y
set style line 1 linecolor rgb "black"
base = "testCase2-cputvscellc-s-"
end = "-T-2000-t-10-m-".muw."-".muo.".data"
plot base."r".end title "Regula Falsi" with linespoints, \
     base."t".end title "Trust Region" with linespoints, \
     base."ta".end title "Trust Region*" with linespoints, \
     base."b".end title "Brent" with linespoints, \
     base."i".end title "Ridder" with linespoints, \
     base."u".end title "Regula Falsi Trust Region" with linespoints ls 1

#pause -1 "Hit some poor key to continue"
