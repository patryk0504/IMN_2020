set terminal png enhanced size 600,300 
set xlabel "x"
set ylabel "y"
set pm3d map
set contour
unset surface
set view map
set cntrparam levels 15
set palette negative

unset key
####################################################################################
#wykresy psi i ksi
#
set title "Psi Q=-1000"
set output "a-1000psi.png" 
splot "fi_c_minus1000.dat" u 1:2:($3!=0?$3:1/0) w l lw 2 palette t ""

set cntrparam levels 30
set title "Ksi Q=-1000"
set output "a-1000ksi.png" 
splot "fi_c_minus1000.dat" u 1:2:($4!=0?$4:1/0) w l lw 2 palette t ""

set title "Psi Q=-4000"
set output "b-4000psi.png" 
splot "fi_c_minus4000.dat" u 1:2:($3!=0?$3:1/0) w l lw 2 palette t ""

set cntrparam levels 30
set title "Ksi Q=-4000"
set output "b-4000ksi.png" 
splot "fi_c_minus4000.dat" u 1:2:($4!=0?$4:1/0) w l lw 2 palette t ""


set title "Psi Q=4000"
set output "c4000psi.png" 
splot "fi_c_plus4000.dat" u 1:2:($3!=0?$3:1/0) w l lw 2 palette t ""

set cntrparam levels 30
set title "Ksi Q=4000"
set output "c4000ksi.png" 
splot "fi_c_plus4000.dat" u 1:2:($4!=0?$4:1/0) w l lw 2 palette t ""

####################################################################################
#wykresy u i v
set title "u Q=-1000"
set output "a-1000ua.png" 
splot "u_v_minus1000.dat" u 1:2:($3!=0?$3:1/0) w pm3d t ""

set title "v Q=-1000"
set output "a-1000va.png" 
splot "u_v_minus1000.dat" u 1:2:($4!=0?$4:1/0) w pm3d t ""

set title "u Q=-4000"
set output "b-4000ua.png" 
splot "u_v_minus4000.dat" u 1:2:($3!=0?$3:1/0) w pm3d t ""

set title "v Q=-4000"
set output "b-4000va.png" 
splot "u_v_minus4000.dat" u 1:2:($4!=0?$4:1/0) w pm3d t ""

set title "u Q=4000"
set output "c4000ua.png" 
splot "u_v_plus4000.dat" u 1:2:($3!=0?$3:1/0) w pm3d t ""

set title "v Q=4000"
set output "c4000va.png" 
splot "u_v_plus4000.dat" u 1:2:($4!=0?$4:1/0) w pm3d t ""