#-----------------------------------------------------------------------------#
#         This GNUPLOT script generates plot the residuals of D and T         #
#-----------------------------------------------------------------------------#
# To use this script, execute : offbeat > log.offbeat & gnuplot Residuals.gp  #
#-----------------------------------------------------------------------------#
# --- Wait until log.offbeat is created :
pause 0.1
# --- Plot setting :
set title "Convergence process"
set xlabel "Iterations"
set ylabel "Residuals"
set logscale y
set terminal push
set terminal pngcairo enhanced font "Times New Roman,14.0" size 900,600
set output 'residuals.png'
# --- Plot residuals for live monitoring :
    
plot "< cat log.offbeat | grep 'Solving for PuFormFactor'  | cut -d' ' -f9 | tr -d ','" title 'Pu' with lines, \
"< cat log.offbeat | grep 'Solving for AmFormFactor'  | cut -d' ' -f9 | tr -d ','" title 'Am' with lines 

set terminal pop
set output 

while (1) {

    # --- Exit if the log.offbeat file it does not update anymore :
    if ( system("cat log.offbeat | grep -v 'end time' | grep -c ' end '") == 1 ) {
        replot
        print "end of simulation"
        pause 3
        quit
    }
    
    # --- Otherwise, re-read:
    unset warnings
    replot
    set warnings
    pause 1

}
