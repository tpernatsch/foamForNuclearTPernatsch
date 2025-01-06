
set y2tics
set ytics nomirror

set xlabel "Time [s]"
set ylabel "Power [W]"
set y2label "External Source Power [W]"

plot file using 1:2 w l axis x1y1 title "Fission",\
    file using 1:3 w l axis x1y1 title "Decay",\
    file using 1:4 w l axis x1y2 title "Ext Src"
