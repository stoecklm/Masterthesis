# ScaFES
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.
#
# Using the following parameters, the corresponding strong scalability diagram
# will be plotted.
# Call: gnuplot <nameOfThisFile>

nameSystem = "ZIH cluster Taurus"
strXlabel = "Number of MPI processes"
sizeFont  = 18
nameBase  = "ZIH_TAURUS-scaling_strong";
strYlabel = "Speedup";
strTitle  = "Strong Scalability"; 
xMin      = 0; 
xMax      = 1024; 
xSteps    = 128; 
yMin      = 0; 
yMax      = 1024; 
ySteps    = 128; 

plotStrongScaling=1
load "scafes-scalingtests.plt";
