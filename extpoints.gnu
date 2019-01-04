set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./EPS.eps"
set title "Array length vs. Search time"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./RESULTS/EPS100' using 1:3 title "EPS 10^2" with linespoints lt 0 pt 5 lw 0.5, \
	'./RESULTS/EPS1000' using 1:3 title "EPS 10^3" with linespoints lt 1 pt 6 lw 0.5, \
	'./RESULTS/EPS10000' using 1:3 title "EPS 10^4" with linespoints lt 2 pt 7 lw 0.5, \
	'./RESULTS/EPS100000' using 1:3 title "EPS 10^5" with linespoints lt 3 pt 8 lw 0.5, \
	'./RESULTS/EPS1000000' using 1:3 title "EPS 10^6" with linespoints lt 4 pt 9 lw 0.5

