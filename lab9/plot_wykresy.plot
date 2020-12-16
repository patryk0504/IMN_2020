set term png


set output "mapaT_100.png"
set xlabel "x"
set ylabel "y"
set title "T, it = 100"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaT_100.dat" i 0 u 1:2:3

set output "mapaT_200.png"
set xlabel "x"
set ylabel "y"
set title "T, it = 200"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaT_200.dat" i 0 u 1:2:3

set output "mapaT_500.png"
set xlabel "x"
set ylabel "y"
set title "T, it = 500"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaT_500.dat" i 0 u 1:2:3

set output "mapaT_1000.png"
set xlabel "x"
set ylabel "y"
set title "T, it = 1000"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaT_1000.dat" i 0 u 1:2:3

set output "mapaT_2000.png"
set xlabel "x"
set ylabel "y"
set title "T, it = 2000"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaT_2000.dat" i 0 u 1:2:3



#######################

set output "mapaDT_100.png"
set xlabel "x"
set ylabel "y"
set title "laplasjan T, it = 100"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaDT_100.dat" i 0 u 1:2:3

set output "mapaDT_200.png"
set xlabel "x"
set ylabel "y"
set title "laplasjan T, it = 200"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaDT_200.dat" i 0 u 1:2:3

set output "mapaDT_500.png"
set xlabel "x"
set ylabel "y"
set title "laplasjan T, it = 500"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaDT_500.dat" i 0 u 1:2:3

set output "mapaDT_1000.png"
set xlabel "x"
set ylabel "y"
set title "laplasjan T, it = 1000"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaDT_1000.dat" i 0 u 1:2:3


set output "mapaDT_2000.png"
set xlabel "x"
set ylabel "y"
set title "laplasjan T, it = 2000"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:40][0:40] "mapaDT_2000.dat" i 0 u 1:2:3