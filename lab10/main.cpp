#include <iostream>
#include <cmath>
#include <stdio.h>

const int nx = 150;
const int nt = 1000;
const double delta = 0.1;
const double dt = 0.05;
const double xa = 7.5;
const double sigma = 0.5;
const double xf = 2.5;

auto *u = new double [nx+1];
auto *u0 = new double [nx+1];
auto *v = new double [nx+1];
auto *vp = new double [nx+1];
auto *a = new double [nx+1];

double af(double x, double t){
    return x == xf ? std::cos(50.*t/(dt*nt)) : 0.;
}

void algorytmVerleta(double alpha, double beta, FILE *fE, FILE *fU,bool isBrzegoweAndPoczatkowe ,bool isInicjalizacjaA){
    if(isBrzegoweAndPoczatkowe){
        //warunki brzegowe
        u[0] = 0;
        u[nx] = 0;
        v[0] = 0;
        v[nx] = 0;

        //warunki poczatkowe
        for(int i=1; i<nx; i++){
            double x = delta*i;
            u[i] = std::exp(-1.0*std::pow(x-xa,2)/(2*sigma*sigma));
            v[i] = 0;
        }
    }
    //zachowanie poprzedniego wyniku
    for(int i=0; i < nx; i++)
        u0[i] = u[i];

    //wyznaczamy a
    if(isInicjalizacjaA) {
        for (int i = 1; i < nx; i++) {
            double el1 = (u[i + 1] - 2. * u[i] + u[i - 1]) / (delta * delta);
            double el2 = beta * (u[i] - u0[i]) / dt;
            double el3 = alpha * af(delta * i, 0);
            a[i] = el1 - el2 + el3;
        }
    }
    //glowna petla algorytmu
    for(int n=1; n<=nt; n++){
        for(int i=0; i<nx; i++){
            vp[i] = v[i] + (dt/2.) * a[i];
            u0[i] = u[i];
        }
        for(int i=0; i<nx; i++){
            u[i] = u[i] + dt*vp[i];
        }
        for(int i=1; i<nx; i++){
            double el1 = (u[i+1] - 2.*u[i] + u[i-1])/(delta*delta);
            double el2 = beta * (u[i] - u0[i])/dt;
            double el3 = alpha * af(delta*i, dt*n);
            a[i] = el1 - el2 + el3;
        }
        for(int i=0; i<nx; i++){
            v[i] = vp[i] + (dt/2.) * a[i];
        }
        double du_dx = 0.;
        double E = 0.;
        for(int i=1; i< nx; i++)
            du_dx += std::pow(v[i],2) + std::pow((u[i+1] - u[i-1])/(2.*delta),2);

        E = 0.25*delta * (std::pow((u[1]-u[0])/delta,2) + std::pow((u[nx]-u[nx-1])/delta,2)) + (0.5*delta) * du_dx;

        //zapis E, u
        fprintf(fE, "%f\t%f\n", dt * n, E);

        for(int i = 0; i < nx; i++)
        {
            fprintf(fU, "%f\t%f\t%f\n", dt*n, delta*i, u[i]);
        }
        fprintf(fU, "\n");
    }
}



int main() {

    //zad1 alpha = 0, beta = 0
    FILE *fE = fopen("E_0_0.dat","w");
    FILE *fU = fopen("U_0_0.dat","w");
    algorytmVerleta(0,0,fE,fU, true, true);
    fclose(fU);
    fclose(fE);
    //alpha = 0, beta 0.1
    fE = fopen("E_0_01.dat","w");
    fU = fopen("U_0_01.dat","w");
    algorytmVerleta(0,0.1,fE,fU, true, true);
    fclose(fU);
    fclose(fE);
    //alpha = 0, beta 1
    fE = fopen("E_0_1.dat","w");
    fU = fopen("U_0_1.dat","w");
    algorytmVerleta(0,1,fE,fU, true, true);
    fclose(fU);
    fclose(fE);

    //zad2 alpha = 1, beta 1
    fE = fopen("E_1_1.dat","w");
    fU = fopen("U_1_1.dat","w");
    algorytmVerleta(1,1,fE,fU, false, false);
    fclose(fU);
    fclose(fE);

    //dealokacja
    delete [] u;
    delete [] u0;
    delete [] v;
    delete [] vp;
    delete [] a;

    return 0;
}
