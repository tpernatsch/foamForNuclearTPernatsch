set terminal png
set output 'sigma.png'
set title "Stress Components"
set xlabel "x (cm)"
set ylabel "Stress (MPa)"
set yrange [-1.2:0.2]
plot "postProcessing/xyGraph/1/radialProfile_sigma_epsilon.xy" using (100*$1):($2/1e6) title "x" with linespoints pointnumber 20, \
     "postProcessing/xyGraph/1/radialProfile_sigma_epsilon.xy" using (100*$1):($5/1e6) title "y" with linespoints pointnumber 21, \
     "postProcessing/xyGraph/1/radialProfile_sigma_epsilon.xy" using (100*$1):($3/1e6) title "xy" with linespoints pointnumber 22

reset
set output 'epsilon.png'
set title "Strain Components"
set xlabel "x (cm)"
set ylabel "Strain (μm/m)"
plot "postProcessing/xyGraph/1/radialProfile_sigma_epsilon.xy" using (100*$1):(1e6*$8) title "x" with linespoints pointnumber 20, \
     "postProcessing/xyGraph/1/radialProfile_sigma_epsilon.xy" using (100*$1):(1e6*$11) title "y" with linespoints pointnumber 21, \
     "postProcessing/xyGraph/1/radialProfile_sigma_epsilon.xy" using (100*$1):(1e6*$9) title "xy" with linespoints pointnumber 22
