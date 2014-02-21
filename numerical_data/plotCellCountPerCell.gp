#set terminal latex
#set terminal wxt
set terminal epslatex size 15cm,8cm color colortext
set output "testCase1-cputvscellcpcell-s-x-T-2000-t-10-viscrat-".muw."_".muo.".tex"
#set format xy "$%g$"
#set format xy "%g"
#set title "CPU time vs time step length"
set xlabel "Cell count"
set ylabel "CPU time [seconds]"
set ylabel rotate left
#set logscale y
#set yrange [ 0.0012  : 0.0009  ]
set auto fix
set style line 1 linecolor rgb "black"
plot "testCase1-cputvscellc-s-r-T-2000-t-10.data" using 1:($2/$1) title "Regula Falsi" with linespoints, \
     "testCase1-cputvscellc-s-u-T-2000-t-10.data" using 1:($2/$1) title "Regula Falsi TR" with linespoints, \
     "testCase1-cputvscellc-s-b-T-2000-t-10.data" using 1:($2/$1) title "Brent" with linespoints, \
     "testCase1-cputvscellc-s-i-T-2000-t-10.data" using 1:($2/$1) title "Ridder" with linespoints, \
     "testCase1-cputvscellc-s-t-T-2000-t-10.data" using 1:($2/$1) title "Trust Region" with linespoints, \
     "testCase1-cputvscellc-s-ta-T-2000-t-10.data" using 1:($2/$1) title "Trust Region*" with linespoints ls 1


#pause -1 "Hit some poor key to continue"
