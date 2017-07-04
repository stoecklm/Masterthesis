# ScaFES
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.
#
# Using the following parameters, the corresponding weak scalability diagram
# will be plotted.
# Call: gnuplot <nameOfThisFile>

nameSystem = "ZIH cluster Taurus"
strXlabel = "Number of MPI processes"
sizeFont  = 18
nameBase  = "ZIH_TAURUS-scaling_weak";
strYlabel = "CPU times in s"; 
strTitle  = "Weak Scalability"; 
xMin      = 0; 
xMax      = 1024; 
xSteps    = 128;
yMin      = 0; 
yMax      = 2; 
ySteps    = 0.25; 

plotStrongScaling=0
load "scafes-scalingtests.plt";
