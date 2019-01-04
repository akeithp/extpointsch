set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./EPSvsEPSP.eps"
set title "Array length vs. Search time"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./RESULTS/EPS' using 1:3 title "EPS" with linespoints lt 1 pt 2 lw 0.5, \
	'./RESULTSP/EPSParallel1' using 1:3 title "EPSP 1 threads" with linespoints lt 2 pt 3 lw 0.5, \
	'./RESULTSP/EPSParallel2' using 1:3 title "EPSP 2 threads" with linespoints lt 3 pt 4 lw 0.5, \
	'./RESULTSP/EPSParallel4' using 1:3 title "EPSP 4 threads" with linespoints lt 4 pt 5 lw 0.5

