# --- Plot setting :
set title "Convergence process"
set xlabel "Iterations"
set ylabel "Residuals"
set logscale y
set terminal push
# --- Plot residuals for live monitoring :
    
plot "< cat log.offbeat | grep 'Solving for T'  | cut -d' ' -f9 | tr -d ','" title 'T' with lines 
set terminal pop
replot

reread