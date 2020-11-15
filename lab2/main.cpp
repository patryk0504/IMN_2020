#include <iostream>
#include <cmath>
#include <stdio.h>
#include <array>

//dane z zadania
const double beta = 0.001;
const int N = 500;
const double y = 0.1;
const int tmax = 100;
const double dt = 0.1;
const double TOL = std::pow(10,-6);
const double alfa = beta * N - y;
double stala_gamma = 0.1;

//funkcje pomocnicze
double funPomoc(double u);
double funF(double u);
///////////////////////////
//zad1
void iteracjePicarda();
//zad2
void iteracjeNewtona();
//zad3
void RK2();

//////////////////////////////////////////////
int main() {
    iteracjePicarda();
    iteracjeNewtona();
    RK2();
    return 0;
}
/////////////////////////////////////////////

double funPomoc(double u){return alfa * u - beta * std::pow(u,2);}
double funF(double u){return (beta*N - stala_gamma)*u - beta * pow(u,2);}

void iteracjePicarda(){
    FILE * fp = fopen("metPicarda.txt","w");
    auto *U = new double [1000];
    //warunek poczatkowy
    U[0] = 1.;
//    fprintf(fp,"%12.6g\t%12.6g\t%12.6g\n",0, U[0],N - U[0]);

    for(int i=1; i<1000; i++){
        double u_prev = U[i-1];
        double u_actual = 0;
        double mi = 0.;
        do{
            u_prev = u_actual;
            u_actual = U[i-1] + (dt/2.) * (funPomoc(U[i-1]) + funPomoc(u_prev));
            mi++;
        }while(mi <= 20 && std::abs(u_actual - u_prev) > TOL);
        U[i] = u_actual;
        fprintf(fp,"%12.6g\t%12.6g\t%12.6g\n",i*dt, U[i],N - U[i]);
    }
    delete [] U;
    fclose(fp);
}

void iteracjeNewtona(){
    FILE * fp = fopen("metNewtona.txt","w");
    auto *U = new double [1000];
    //warunek poczatkowy
    U[0] = 1.;
//    fprintf(fp,"%12.6g\t%12.6g\t%12.6g\n",0, U[0],N - U[0]);

    for(int i=1; i<1000; i++){
        double u_prev = U[i-1];
        double u_actual = 0;
        double mi = 0.;
        do{
            u_prev = u_actual;
            u_actual = u_prev - (u_prev - U[i-1] - (dt/2.) * ( funPomoc(U[i-1]) + funPomoc(u_prev) ))/(1-(dt/2.)*(alfa - 2*beta*u_prev));
            mi++;
        }while(mi <= 20 && std::abs(u_actual - u_prev) > TOL);
        U[i] = u_actual;
        fprintf(fp,"%12.6g\t%12.6g\t%12.6g\n",i*dt, U[i],N - U[i]);
    }
    delete [] U;
    fclose(fp);
}

void RK2(){
    FILE * fp = fopen("metRK2.txt","w");
    auto *U = new double [1000];
    auto *resultOfRK2 = new double [1001];

    //warunki poczatkowe
    U[0] = 1.;
    resultOfRK2[0]=1.;
//    fprintf(fp,"%12.6g\t%12.6g\t%12.6g\n",0, U[0],N - U[0]);

    //tablica Butchera
    //maierz A
    std::array<std::array<double,2>,2> a = {{{0.25, 0.25 - sqrt(3)/6},{0.25 + sqrt(3)/6,0.25}}};
    //wektor B (wartosci takie same)
    double b = 0.5;
    //wektor C
    std::array<double,2> c = {{0.5 - sqrt(3)/6, 0.5 + sqrt(3)/6}};
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(int i=1; i<1000; i++) {
        auto *U1 = new double[1001];
        auto *U2 = new double[1001];
        U1[0] = resultOfRK2[i - 1];
        U2[0] = resultOfRK2[i - 1];
        int mi = 0;
        do {
            //wektor F
            double F1 = U1[mi] - resultOfRK2[i - 1] - (dt * (a[0][0] * funPomoc(U1[mi]) + a[0][1] * funPomoc(U2[mi])));
            double F2 = U2[mi] - resultOfRK2[i - 1] - (dt * (a[1][0] * funPomoc(U1[mi]) + a[1][1] * funPomoc(U2[mi])));
            //wspolczynniki macierzy M
            double m11 = 1 - dt * a[0][0] * (alfa - 2 * beta * U1[mi]);
            double m12 = -dt * a[0][1] * (alfa - 2 * beta * U2[mi]);
            double m21 = -dt * a[1][0] * (alfa - 2 * beta * U1[mi]);
            double m22 = 1 - dt * a[1][1] * (alfa - 2 * beta * U2[mi]);
            //wektor delta U
            double dU1 = (F2 * m12 - F1 * m22) / (m11 * m22 - m12 * m21);
            double dU2 = (F1 * m21 - F2 * m11) / (m11 * m22 - m12 * m21);
            //zapisujemy kolejne przyblizenia
            U1[mi + 1] = U1[mi] + dU1;
            U2[mi + 1] = U2[mi] + dU2;
            mi++;
//            std::cout<<mi<<std::endl;
        } while (mi <= 20 && U1[mi] > TOL && U2[mi] > TOL);

//        std::cout<<"index1: "<<index1<<"    index2: "<<index2<<std::endl;

        //rozwiazanie w chwili n+1
        resultOfRK2[i] = resultOfRK2[i - 1] + dt * ((b * funF(U1[mi])) + b * funF(U2[mi]));
        fprintf(fp, "%12.6g\t%12.6g\t%12.6g\n", dt * i, resultOfRK2[i], N - resultOfRK2[i]);

        delete [] U1;
        delete [] U2;
    }
    fclose(fp);
}
