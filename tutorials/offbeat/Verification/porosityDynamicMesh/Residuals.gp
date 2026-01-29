set logscale y
set title "Residuals"
set ylabel 'Residual'
set xlabel 'Iteration'
plot    "< cat log.offbeat | grep 'Solving for T'  | cut -d' ' -f9 | tr -d ','" title 'T' with lines,\
        "< cat log.offbeat | grep 'Solving for porosity' | cut -d' ' -f9 | tr -d ','" title 'porosity' with lines
pause 1
reread