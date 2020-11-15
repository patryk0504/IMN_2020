set term png

set logscale x

set title "Relaksacja Globalna"
set out "relaksacja_globalna_S.png"
set xlabel "iteracje"
set ylabel "S"
plot "S_06_global.txt" u 1:2 w l lw 1 t 'w=0.6', "S_1_global.txt" u 1:2 w l lw 1 t 'w=1.0'

set title "Relaksacja Lokalna"
set out "relaksacja_lokalna_S.png"
set xlabel "iteracje"
set ylabel "S"
plot "S_1_lokal.txt" u 1:2 w l lw 1 t 'wt=1.0', "S_1_4_lokal.txt" u 1:2 w l lw 1 t 'wt=1.4', "S_1_8_lokal.txt" u 1:2 w l lw 1 t 'wt=1.8', "S_1_9_lokal.txt" u 1:2 w l lw 1 t 'wt=1.9'

unset logscale x

set xrange [0.0:15.0]
set yrange [0.0:10.0]

set title "Relaksacja Globalna - potencial V (0.6)"
set out "globalna_v_0_6.png"
set xlabel "X"
set ylabel "Y"
plot "V_06_global.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Relaksacja Globalna - potencial V (1.0)"
set out "globalna_v_1_0.png"
set xlabel "X"
set ylabel "Y"
plot "V_1_global.txt" u 1:2:3 w p pt 7 palette t 'V'

set title "Relaksacja globalna - blad (w = 0.6)"
set out "globalna_blad_0_6.png"
set cbrange [ 0.00000 : 0.003 ]
set xlabel "X"
set ylabel "Y"
plot "Err_06_global.txt" u 1:2:3 w p pt 7 palette t 'Blad'

set title "Relaksacja globalna - blad (w = 1.0)"
set out "globalna_blad_1_0.png"
set cbrange [ 0.00000 : 0.0025 ]
set xlabel "X"
set ylabel "Y"
plot "Err_1_global.txt" u 1:2:3 w p pt 7 palette t 'Blad'