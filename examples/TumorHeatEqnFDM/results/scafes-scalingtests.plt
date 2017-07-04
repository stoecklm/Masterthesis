# ScaFES
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.
#
# Plots a weak or a strong scalability diagram depending on the parameter
# 'plotStrongScaling'.

nameDataFile=nameBase.".dat"

if (plotStrongScaling==1) \
    reset; \
    plot nameDataFile using 1:2; \
    yPeak = GPVAL_DATA_Y_MAX

set terminal wxt persist
set object 1 rectangle from screen 0,0 to screen 1,1 \
           fillcolor rgb"#eee8cd" behind
set style line 1 lt rgb "#000080" lw 3
set style line 2 lt 1 lw 1
set format x '%.0f'
set format y '%8.1f'
set xtics xMin,xSteps,xMax
set ytics yMin,ySteps,yMax
set xrange [xMin:xMax]
set yrange [yMin:yMax]
set xlabel strXlabel font 'Helvetica, '.sizeFont
set ylabel strYlabel font 'Helvetica, '.sizeFont
set size square {1,1}
set grid
set title strTitle font 'Helvetica, '.sizeFont

if (plotStrongScaling==1) \
    plot nameDataFile using  ($1):(yPeak/($2)) t nameSystem w lp ls 1; \
else \
    plot nameDataFile using  ($1):($2) t nameSystem w lp ls 1;

set terminal postscript color 'Helvetica, '.sizeFont
set output nameBase.".eps"
replot

set terminal png
set output nameBase.".png"
replot
