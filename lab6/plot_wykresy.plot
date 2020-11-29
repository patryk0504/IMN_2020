set term png

set output "mapa50_50.png"
set xlabel "x"
set ylabel "y"
set title "nx=50 ny=50"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:5][0:5] "map50_50.txt" i 0 u 1:2:3

reset

set output "mapa100_100.png"
set xlabel "x"
set ylabel "y"
set title "nx=100 ny=100"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:10][0:10] "map100_100.txt" i 0 u 1:2:3

reset

set output "mapa200_200.png"
set xlabel "x"
set ylabel "y"
set title "nx=200 ny=200"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:20][0:20] "map200_200.txt" i 0 u 1:2:3

reset

set output "mapa(e1_1_e2_1).png"
set xlabel "x"
set ylabel "y"
set title "nx=ny=100, eps1 = 1, eps2 = 1"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:10][0:10] "map_e_1_1.txt" i 0 u 1:2:3

reset

set output "mapa(e1_1_e2_2).png"
set xlabel "x"
set ylabel "y"
set title "nx=ny=100, eps1 = 1, eps2 = 2"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1
set cbrange [-0.8:0.8]


splot [0:10][0:10] "map_e_1_2.txt" i 0 u 1:2:3

reset

set output "mapa(e1_1_e2_10).png"
set xlabel "x"
set ylabel "y"
set title "nx=ny=100, eps1 = 1, eps2 = 10"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1
set cbrange [-0.8:0.8]


splot [0:10][0:10] "map_e_1_10.txt" i 0 u 1:2:3