#set terminal latex
#set terminal wxt
set terminal epslatex size 15cm,8cm color colortext
set output "testCase1-cputvscellc-s-x-T-2000-t-10.tex"
set format xy "$%g$"
#set format xy "%g"
#set title "CPU time vs time step length"
set xlabel "Cell count [days]"
set ylabel "CPU time [seconds]"
set ylabel rotate left
set logscale y
plot "testCase1-cputvscellc-s-r-T-2000-t-10.data" title "Regula Falsi" with linespoints, \
     "testCase1-cputvscellc-s-b-T-2000-t-10.data" title "Brent" with linespoints, \
     "testCase1-cputvscellc-s-i-T-2000-t-10.data" title "Ridder" with linespoints, \
     "testCase1-cputvscellc-s-t-T-2000-t-10.data" title "Trust Region" with linespoints

pause -1 "Hit some poor key to continue"
