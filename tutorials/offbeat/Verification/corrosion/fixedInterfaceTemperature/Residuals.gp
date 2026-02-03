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
    
plot "< cat log.offbeat | grep 'Solving for Dx' | cut -d' ' -f9 | tr -d ','" title 'Dx' with lines, \
"< cat log.offbeat | grep 'Solving for Dy' | cut -d' ' -f9 | tr -d ','" title 'Dy' with lines, \
"< cat log.offbeat | grep 'Solving for Dz' | cut -d' ' -f9 | tr -d ','" title 'Dz' with lines, \
"< cat log.offbeat | grep 'Solving for DDx' | cut -d' ' -f9 | tr -d ','" title 'DDx' with lines, \
"< cat log.offbeat | grep 'Solving for DDy' | cut -d' ' -f9 | tr -d ','" title 'DDy' with lines, \
"< cat log.offbeat | grep 'Solving for DDz' | cut -d' ' -f9 | tr -d ','" title 'DDz' with lines, \
"< cat log.offbeat | grep 'Solving for T'  | cut -d' ' -f9 | tr -d ','" title 'T' with lines
set terminal pop
set output 
replot
# --- Exit if the log.offbeat file it does not update anymore :
if ( ` cat log.offbeat | grep -c " end " ` == 1){
    pause 1
    quit
}
pause 1
# --- Otherwise, re-read:
reread