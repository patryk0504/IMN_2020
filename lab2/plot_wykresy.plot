set term png
set xlabel "Iteracje"
set ylabel "Wartosci"

set title "Iteracja Picarda"
set out "IterPicard.png"
plot "metPicarda.txt" u 1:2 w l lw 3 t 'u(t)',"metPicarda.txt" u 1:3 w l lw 2 t 'N - u(t)'

set title "Iteracja Newtona"
set out "IterNewton.png"
plot "metNewtona.txt" u 1:2 w l lw 3 t 'u(t)',"metNewtona.txt" u 1:3 w l lw 2 t 'N - u(t)'

set title "Metoda RK2"
set out "RK2.png"
plot "metRK2.txt" u 1:2 w l lw 3 t 'u(t)',"metRK2.txt" u 1:3 w l lw 2 t 'N - u(t)'