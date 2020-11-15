set term jpeg

set xlabel "iter"
set ylabel "S(iter)"

# set logscale y
set yrange [4.2:5.6]

set title "Relaksacja wielosiatkowa"
set out "S.jpeg"
plot "S_16.txt" u 1:2 w l lw 1 t 'k=16',"S_8.txt" u 1:2 w l lw 1 t 'k=8',"S_4.txt" u 1:2 w l lw 1 t 'k=4',"S_2.txt" u 1:2 w l lw 1 t 'k=2',"S_1.txt" u 1:2 w l lw 1 t 'k=1'

# unset logscale y
unset yrange

set xlabel "x"
set ylabel "y"
set zlabel "w"

set logscale z

set title "Relaksacja wielosiatkowa"
set out "V_k_16.jpeg"
plot [-3:29][-3:29] "V_16.txt" u 1:2:3 w p pt 5 ps 8 palette notitle

set title "Relaksacja wielosiatkowa"
set out "V_k_8.jpeg"
plot [-3:29][-3:29] "V_8.txt" u 1:2:3 w p pt 5 ps 4 palette notitle

set title "Relaksacja wielosiatkowa"
set out "V_k_4.jpeg"
plot [-3:29][-3:29] "V_4.txt" u 1:2:3 w p pt 5 ps 2 palette notitle

set title "Relaksacja wielosiatkowa"
set out "V_k_2.jpeg"
plot [-3:29][-3:29] "V_2.txt" u 1:2:3 w p pt 5 ps 1 palette notitle

set title "Relaksacja wielosiatkowa"
set out "V_k_1.jpeg"
plot [-3:29][-3:29] "V_1.txt" u 1:2:3 w p pt 7 palette notitle