
set yr [300:1500]
set xlabel "Time [s]"
set ylabel "Temperature [K]"
plot file using 1:2 w l title "TFuel",\
    file using 1:3 w l title "TClad",\
    file using 1:4 w l title "TCool",\
