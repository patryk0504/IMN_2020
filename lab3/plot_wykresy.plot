set term png



set title "Metoda Trapezow - x(t)"
set out "met_trapez_xt.png"
plot "trapezy_1_TOL_5.txt" u 1:3 w l lw 1 t 'tol = 10^{-5}', "trapezy_1_TOL_2.txt" u 1:3 w l lw 1 t 'tol = 10^{-2}'

set title "Metoda Trapezow - v(t)"
set out "met_trapez_vt.png"
plot "trapezy_1_TOL_5.txt" u 1:4 w l lw 1 t 'tol = 10^{-5}', "trapezy_1_TOL_2.txt" u 1:4 w l lw 1 t 'tol = 10^{-2}'

set title "Metoda Trapezow - dt(t)"
set out "met_trapez_dtt.png"
plot "trapezy_1_TOL_5.txt" u 1:2 w l lw 1 t 'tol = 10^{-5}', "trapezy_1_TOL_2.txt" u 1:2 w l lw 1 t 'tol = 10^{-2}'

set title "Metoda Trapezow - v(x)"
set out "met_trapez_vx.png"
plot "trapezy_1_TOL_5.txt" u 3:4 w l lw 1 t 'tol = 10^{-5}', "trapezy_1_TOL_2.txt" u 3:4 w l lw 1 t 'tol = 10^{-2}'


set title "Metoda RK2 - x(t)"
set out "met_rk2_xt.png"
plot "rk2_1_TOL_5.txt" u 1:3 w l lw 1 t 'tol = 10^{-5}', "rk2_1_TOL_2.txt" u 1:3 w l lw 1 t 'tol = 10^{-2}'

set title "Metoda RK2 - v(t)"
set out "met_rk2_vt.png"
plot "rk2_1_TOL_5.txt" u 1:4 w l lw 1 t 'tol = 10^{-5}', "rk2_1_TOL_2.txt" u 1:4 w l lw 1 t 'tol = 10^{-2}'

set title "Metoda RK2 - dt(t)"
set out "met_rk2_dtt.png"
plot "rk2_1_TOL_5.txt" u 1:2 w l lw 1 t 'tol = 10^{-5}', "rk2_1_TOL_2.txt" u 1:2 w l lw 1 t 'tol = 10^{-2}'

set title "Metoda RK2 - v(x)"
set out "met_rk2_vx.png"
plot "rk2_1_TOL_5.txt" u 3:4 w l lw 1 t 'tol = 10^{-5}', "rk2_1_TOL_2.txt" u 3:4 w l lw 1 t 'tol = 10^{-2}'