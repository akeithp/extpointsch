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
	'./RESULTS/EPS1' using 1:3 title "EPS Normal" with linespoints lt 1 pt 2 lw 0.5, \
	'./RESULTS/EPS0' using 1:3 title "EPS Uniform" with linespoints lt 2 pt 3 lw 0.5

