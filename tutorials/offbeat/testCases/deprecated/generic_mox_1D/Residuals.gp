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

# --- Plot residuals for live monitoring :
    
plot "< cat log.offbeat | grep 'Solving for Dx' | cut -d' ' -f9 | tr -d ','" title 'Dx' with lines, \
"< cat log.offbeat | grep 'Solving for Dz' | cut -d' ' -f9 | tr -d ','" title 'Dz' with lines, \
"< cat log.offbeat | grep 'Solving for T'  | cut -d' ' -f9 | tr -d ','" title 'T' with lines, \
"< cat log.offbeat | grep 'Solving for neutronFlux0'  | cut -d' ' -f9 | tr -d ','" title 'neutronFlux0' with lines, \
"< cat log.offbeat | grep 'Solving for PuFormFactor'  | cut -d' ' -f9 | tr -d ','" title 'PuFormFactor' with lines, \
"< cat log.offbeat | grep 'Solving for AmFormFactor'  | cut -d' ' -f9 | tr -d ','" title 'AmFormFactor' with lines, \
"< cat log.offbeat | grep 'Solving for porosity'  | cut -d' ' -f9 | tr -d ','" title 'porosity' with lines

pause 1
# --- Otherwise, re-read:
reread