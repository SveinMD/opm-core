set terminal latex
#set terminal wxt
set terminal epslatex size 15cm,8cm color colortext
set output "testCase1-error-s-x-T-2-t-0_1.tex"
set format xy "$%g$"
#set format xy "%g"
#set title "CPU time vs time step length"
set xlabel "Absolute error"
set ylabel "Time step"
set ylabel rotate left
set logscale y
set style line 1 linecolor rgb "black"
base = "diff_s_"
end = ".data"
plot base."t".end title "Trust Region" with linespoints, \
     base."ta".end title "Trust Region*" with linespoints, \
     base."b".end title "Brent" with linespoints, \
     base."i".end title "Ridder" with linespoints
