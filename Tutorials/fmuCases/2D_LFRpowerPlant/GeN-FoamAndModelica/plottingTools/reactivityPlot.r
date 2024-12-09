
set xlabel "Time [s]"
set ylabel "Reactivity [pcm]"
set yr [-100:100]
plot file using 1:2 w l title "Total",\
    file using 1:4 w l title "TFuel",\
    file using 1:3 w l title "Doppler",\
    file using 1:5 w l title "TClad",\
    file using 1:6 w l title "TCool",\
    file using 1:7 w l title "rhoCool",\
    file using 1:8 w l title "TStruct",\
    file using 1:9 w l title "driveline",\
    file using 1:10 w l title "Boron",\
