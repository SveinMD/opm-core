#set terminal latex
#set terminal wxt
set terminal epslatex size 15cm,8cm color colortext
set output "testCase1-cputvsdt-s-x-T-2000.tex"
set format xy "$%g$"
#set format xy "%g"
set title "CPU time vs time step length"
set xlabel "Time step [days]"
set ylabel "CPU time [seconds]"
set ylabel rotate left
set logscale y
plot "testCase1-cputvsdt-s-r-T-2000.data" title "Regula Falsi" with linespoints, \
     "testCase1-cputvsdt-s-b-T-2000.data" title "Brent" with linespoints, \
     "testCase1-cputvsdt-s-i-T-2000.data" title "Ridder" with linespoints, \
     "testCase1-cputvsdt-s-t-T-2000.data" title "Trust Region" with linespoints

pause -1 "Hit some poor key to continue"
