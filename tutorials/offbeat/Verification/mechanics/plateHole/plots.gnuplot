set terminal png
set output 'sigma.png'
set title "Stress Components"
set xlabel "y (cm)"
set ylabel "Stress (MPa)"
plot [0.005:0.04] [0:] "postProcessing/xyGraph/1/symmetry_sigma_epsilon.xy" using ($1):($2/1e6) title "xx" with linespoints, (1+(0.125/((100*x)**2))+(0.09375/((100*x)**4))) title "analytical"

set output 'sigma_norm.png'
set title "Normalised Stress"
set xlabel "y (cm)"
set ylabel "Normalised Stress (-)"
plot [0.005:0.04] [0:] "postProcessing/xyGraph/1/symmetry_sigma_epsilon.xy" using ($1):(7.0/8.0*$2/1e6) title "xx" with linespoints, 7.0/8.0 * (1+(0.125/((100*x)**2))+(0.09375/((100*x)**4))) title "analytical"

set output 'epsilon.png'
set title "Strain Components"
set xlabel "y (cm)"
set ylabel "Strain x 10^{-6} (-)"
plot "postProcessing/xyGraph/1/symmetry_sigma_epsilon.xy" using ($1):($8*1e6) title "xx" with lines
