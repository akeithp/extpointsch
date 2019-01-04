set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./extpointsP.eps"
set title "Array length vs. Seach time Parallel"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./resultsP/extpoints1BLK1000' using 1:4 title "EPSP1 BLK1000" with linespoints lt 0 pt 7 lw 0.5, \
	'./resultsP/extpoints2BLK1000' using 1:4 title "EPSP2 BLK1000" with linespoints lt 1 pt 8 lw 0.5,\
	'./resultsP/extpoints4BLK1000' using 1:4 title "EPSP4 BLK1000" with linespoints lt 2 pt 9 lw 0.5,\
	'./resultsP/extpoints1BLK20000' using 1:4 title "EPSP1 BLK20000" with linespoints lt 3 pt 5 lw 0.5, \
	'./resultsP/extpoints2BLK20000' using 1:4 title "EPSP2 BLK20000" with linespoints lt 4 pt 6 lw 0.5, \
	'./resultsP/extpoints4BLK20000' using 1:4 title "EPSP4 BLK20000" with linespoints lt 5 pt 4 lw 0.5

