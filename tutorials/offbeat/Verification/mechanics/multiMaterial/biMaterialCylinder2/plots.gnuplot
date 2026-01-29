set terminal png
set output 'sigma.png'
set title "Stress Components"
set xlabel "Radius (cm)"
set ylabel "Stress (MPa)"
plot "postProcessing/xyGraph/1/radialProfile_sigmaCyl_epsilon.xy" using (100*$1):($2/1e6) title "Radial" with linespoints pointnumber 20, \
     "postProcessing/xyGraph/1/radialProfile_sigmaCyl_epsilon.xy" using (100*$1):($5/1e6) title "Hoop" with linespoints pointnumber 21, \
     "postProcessing/xyGraph/1/radialProfile_sigmaCyl_epsilon.xy" using (100*$1):($3/1e6) title "Shear" with linespoints pointnumber 22

set output 'epsilon.png'
set title "Strain Components"
set xlabel "Radius (cm)"
set ylabel "Strain (-)"
plot "postProcessing/xyGraph/1/radialProfile_sigmaCyl_epsilon.xy" using (100*$1):8 title "Radial" with linespoints pointnumber 20, \
     "postProcessing/xyGraph/1/radialProfile_sigmaCyl_epsilon.xy" using (100*$1):11 title "Hoop" with linespoints pointnumber 21, \
     "postProcessing/xyGraph/1/radialProfile_sigmaCyl_epsilon.xy" using (100*$1):9 title "Shear" with linespoints pointnumber 22

set output 'D.png'
set title "Displacement"
set xlabel "Radius (cm)"
set ylabel "Displacement (um)"
plot "postProcessing/xyGraph/1/radialProfile_D.xy" using (100*$1):(1e6*$2) title "Radial" with linespoints pointnumber 20, \
     "postProcessing/xyGraph/1/radialProfile_D.xy" using (100*$1):(1e6*$3) title "Azimuthal" with linespoints pointnumber 21
