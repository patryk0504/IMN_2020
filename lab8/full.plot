set term png enhanced size 800,600 

set size ratio -1

set output "vx_map.png"
set view map
set title "vx"
splot "resVx.dat" using 1:2:3 with pm3d

set output "vy_map.png"
set view map
set title "vy"
splot "resVy.dat" using 1:2:3 with pm3d

set terminal png size 800,600
set size ratio 1

set xlabel "tn"
set ylabel "C(tn)"

set yrange [0.5:1.05]
set title "Calka gestosci - c(tn)"
set out "c(tn).png"
plot "res.dat" u 1:2 w l lw 1 t 'C(tn) : D = 0', "res2.dat" u 1:2 w l lw 1 t "C(tn) : D = 0.1"

set xlabel "tn"
set ylabel "Xsr(tn)"

set yrange [0:4]
set title "Srednie polozenie - xsr(tn)"
set out "xsr(tn).png"
plot "res.dat" u 1:3 w l lw 1 t 'Xsr(tn) : D = 0', "res2.dat" u 1:3 w l lw 1 t 'Xsr(tn) : D = 0.1'


reset
set term gif size 800,400 animate delay 10
set output "anim1_D_0.gif"
n=30    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]

do for [i=0:n] {
  file = sprintf("mapa_%i.dat",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 

reset
set term gif size 800,400 animate delay 10
set output "anim2_D_0_1.gif"
n=30    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]
do for [i=0:n] {
  file = sprintf("mapa2_%i.dat",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 