set term png
set xlabel "Iteracje"
set ylabel "Wartosci"

set title "Euler"
set out "euler.png"
plot "euler_1.txt" u 1:2 w l lw 3 t 'dt = 0.01', "euler_2.txt" u 1:2 w l lw 2 t 'dt = 0.1', "euler_3.txt" u 1:2 w l lw 2 t 'dt = 1.0', "metAnalityczna.txt" u 1:2 w p pt 2 ps 0.2 t 'dok'

set title "RK2"
set out "RK2.png"
plot "rk2_1.txt" u 1:2 w l lw 4 t 'dt = 0.01', "rk2_2.txt" u 1:2 w l lw 2 t 'dt = 0.1', "rk2_3.txt" u 1:2 w l lw 2 t 'dt = 1.0', "metAnalityczna.txt" u 1:2 w p pt 2 ps 0.2 t 'dok'

set title "RK4"
set out "RK4.png"
plot "rk4_1.txt" u 1:2 w l lw 4 t 'dt = 0.01', "rk4_2.txt" u 1:2 w l lw 2 t 'dt = 0.1', "rk4_3.txt" u 1:2 w l lw 2 t 'dt = 1.0', "metAnalityczna.txt" u 1:2 w p pt 2 ps 0.2 t 'dok'


set title "RLC omega"
set out "rlc_omega.png"
plot "rlc_0_5.txt" u 1:2 w l lw 1 t 'Wv = 0.5*w0', "rlc_0_8.txt" u 1:2 w l lw 2 t 'Wv = 0.8*w0', "rlc_1.txt" u 1:2 w l lw 3 t 'Wv = dt = 1.0 * w0', "rlc_1_2.txt" u 1:2 w l lw 4 t "Wv = 1.2 * w0"

set title "RLC I"
set out "rlc_I.png"
plot "rlc_0_5.txt" u 1:3 w l lw 1 t 'Wv = 0.5*w0', "rlc_0_8.txt" u 1:3 w l lw 2 t 'Wv = 0.8*w0', "rlc_1.txt" u 1:3 w l lw 3 t 'Wv = dt = 1.0 * w0', "rlc_1_2.txt" u 1:3 w l lw 4 t "Wv = 1.2 * w0"

############################################## błędy
# set logscale y

set title "Euler - blad"
set out "euler_blad.png"
plot "euler_1_blad.txt" u 1:2 w l lw 4 t 'dt = 0.01', "euler_2_blad.txt" u 1:2 w l lw 4 t 'dt = 0.1', "euler_3_blad.txt" u 1:2 w l lw 3 t 'dt = 1.0'

set title "RK2 - blad"
set out "RK2_blad.png"
plot "rk2_1_blad.txt" u 1:2 w l lw 4 t 'dt = 0.01', "rk2_2_blad.txt" u 1:2 w l lw 4 t 'dt = 0.1', "rk2_3_blad.txt" u 1:2 w l lw 3 t 'dt = 1.0'

set title "RK4 - blad"
set out "RK4_blad.png"
plot "rk4_1_blad.txt" u 1:2 w l lw 4 t 'dt = 0.01', "rk4_2_blad.txt" u 1:2 w l lw 4 t 'dt = 0.1', "rk4_3_blad.txt" u 1:2 w l lw 3 t 'dt = 1.0'