set term png



set xlabel "t"
set ylabel "E"
set view map


set out 'E_0_0.png'
set title "E(t) alpha=0 beta=0"
plot 'E_0_0.dat' u 1:2 w l t ""

set out 'E_0_01.png'
set title "E(t) alpha=0 beta=0.1"
plot 'E_0_01.dat' u 1:2 w l t ""

set out 'E_0_1.png'
set title "E(t) alpha=0 beta=1"
plot 'E_0_1.dat' u 1:2 w l t ""

set out 'E_1_1.png'
set title "E(t) alpha=1 beta=1"
plot 'E_1_1.dat' u 1:2 w l t ""

set xlabel "t"
set ylabel "x"
set view map

set out 'mapa_0_0.png'
set title "alpha=0 beta=0"
splot 'U_0_0.dat' u 1:2:3 w pm3d t ""

set out 'mapa_0_01.png'
set title "alpha=0 beta=0.1"
splot 'U_0_01.dat' u 1:2:3 w pm3d t ""

set out 'mapa_0_1.png'
set title "alpha=0 beta=1"
splot 'U_0_1.dat' u 1:2:3 w pm3d t ""

set out 'mapa_1_1.png'
set title "alpha=1 beta=1"
splot 'U_1_1.dat' u 1:2:3 w pm3d t ""