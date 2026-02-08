set terminal png
set output 'sigma_xx.png'
set title "Bending Stress"
set xlabel "y (m)"
set ylabel "Stress (kPa)"
plot [0:] [:] "postProcessing/xyGraph/1/line1_sigma.xy" using ($1):($2/1e3) title "x=0.8" with lines, \
              "postProcessing/xyGraph/1/line2_sigma.xy" using ($1):($2/1e3) title "x=2.4" with lines, \
              "postProcessing/xyGraph/1/line3_sigma.xy" using ($1):($2/1e3) title "x=3.8" with lines, \
              [:0.54] 438.96*(1 - 2*x/0.54) title "analytical", \
              [:0.42] 362.81*(1 - 2*x/0.42) title "analytical", \
              [:0.315] 80.62*(1 - 2*x/0.315) title "analytical"

set output 'sigma_xy.png'
set title "Shear Stress"
set xlabel "y (m)"
set ylabel "Stress (kPa)"
plot "postProcessing/xyGraph/1/line1_sigma.xy" using ($1):($3/1e3)  title "x=0.8" with lines, \
     "postProcessing/xyGraph/1/line2_sigma.xy" using ($1):($3/1e3)  title "x=2.4" with lines, \
     "postProcessing/xyGraph/1/line3_sigma.xy" using ($1):($3/1e3)  title "x=3.8" with lines, \
     [:0.54] -60.9663*(3*(x**2)/0.54 - 2*x) + 254.0263*((x**2) - 0.54*x)  	title "analytical", \
	 [:0.42] -64.7878*(3*(x**2)/0.42 - 2*x) + 539.8985*((x**2) - 0.42*x)  title "analytical", \
	 [:0.315] -19.1964*(3*(x**2)/0.315 - 2*x) + 1279.7594*((x**2) - 0.315*x) title "analytical"
