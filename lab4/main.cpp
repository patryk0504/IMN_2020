#include <iostream>
#include <cmath>
#include <stdio.h>
#include <array>
#include <thread>

//przyjete wartosci parametrow w zadaniu
double eps = 1.;
double delta = 0.1;
const int nx = 150;
const int ny = 100;

double V1 = 10.;
double V2 = 0.;

double xmax = delta*nx;
double ymax = delta*ny;

double sigmax = 0.1 * xmax;
double sigmay = 0.1 * ymax;

double TOL = std::pow(10,-8);

//funkcje gestosci
double gestosc1(double x, double y){
    double wykladnik = -1.*std::pow((x-0.35*xmax),2)/(sigmax*sigmax) - 1.*std::pow((y-0.5*ymax),2)/(sigmay*sigmay);
    return std::exp(wykladnik);
}

double gestosc2(double x, double y){
    double wykladnik = -1.*std::pow((x-0.65*xmax),2)/(sigmax*sigmax) - 1.*std::pow((y-0.5*ymax),2)/(sigmay*sigmay);
    return -std::exp(wykladnik);
}

//metody
void globalna(double wg, std::string filename_S, std::string filename_V, std::string filename_Err);
void lokalna(double wl, std::string filename_S);

int main() {

    //zwykle wywolania (jezeli bylby problem z std::thread)
    std::cout<<"Program is running... Please wait.\n"<<std::endl;

   globalna(0.6, "S_06_global.txt", "V_06_global.txt", "Err_06_global.txt");
   globalna(1.0, "S_1_global.txt", "V_1_global.txt", "Err_1_global.txt");
   lokalna(1.0,"S_1_lokal.txt");
   lokalna(1.4,"S_1_4_lokal.txt");
   lokalna(1.8,"S_1_8_lokal.txt");
   lokalna(1.9,"S_1_9_lokal.txt");


////std::thread
    // std::cout<<"Program is running... Please wait.\n"<<std::endl;

    // std::thread th1 (globalna,0.6, "S_06_global.txt", "V_06_global.txt", "Err_06_global.txt");
    // std::thread th2 (globalna,1.0, "S_1_global.txt", "V_1_global.txt", "Err_1_global.txt");
    // std::thread th3 (lokalna,1.0,"S_1_lokal.txt" );
    // std::thread th4 (lokalna,1.4,"S_1_4_lokal.txt");
    // th1.join();
    // th2.join();
    // th3.join();
    // th4.join();
    // std::thread th5 (lokalna,1.8,"S_1_8_lokal.txt" );
    // std::thread th6 (lokalna,1.9,"S_1_9_lokal.txt" );
    // th5.join();
    // th6.join();
//////////////////////////////////////////////////////////////////////////////////////////////////
    return 0;
}


void globalna(double wg, std::string filename_S, std::string filename_V, std::string filename_Err){

    std::cout<<"Relaksacja globalna dla wg = " +  std::to_string(wg) + " ...\n";

    FILE *fp_S = fopen(filename_S.c_str(),"w");
    FILE *fp_V = fopen(filename_V.c_str(),"w");
    FILE *fp_Err = fopen(filename_Err.c_str(),"w");

    double Vn [nx+5][nx+5] = {0.};
    double P [nx+5][ny+5] = {0.};
    double Vs [nx+5][ny+5] = {0.};

    //tablica gestosci i wstepne rozwiazania Vs Vn
    for(int i=0; i<nx+1; i++){
        Vs[i][0] = V1;
        Vn[i][0] = V1;
        for(int j=0; j<ny+1; j++){
            P[i][j] = gestosc1(i*delta, j*delta) + gestosc2(i*delta, j*delta);
        }
    }

    int count = 0;
    double S[100000] = {0.};
    S[0] = 0.;

    do{
        count++;
        //I etap - elementy poza brzegowymi
        for(int i=1; i<nx; i++){
            for(int j = 1; j<ny; j++){
                Vn[i][j] = 1/4. * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + (delta*delta)/eps * P[i][j]);
//                std::cout<<Vn[i][j]<<std::endl;
            }
        }
        //II etap - WB Neumanna
        for(int j=1; j<ny+1; j++){
            Vn[0][j] = Vn[1][j];
            Vn[nx][j] = Vn[nx-1][j];
        }
        //III etap - mieszamy rozwiazania
        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny+1; j++){
                Vs[i][j]=(1.-wg) * Vs[i][j] + wg * Vn[i][j];
            }
        }
        //warunek stopu
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                S[count] += (delta*delta) * (0.5 * std::pow((Vn[i+1][j] - Vn[i][j])/delta, 2) + 0.5*std::pow((Vn[i][j+1]-Vn[i][j])/delta, 2 ) - (P[i][j]*Vn[i][j]) );
            }
        }
        fprintf(fp_S,"%d\t%12.6g\n",count-1,S[count]);
    }while(std::abs((S[count] - S[count-1]) / S[count-1]) > TOL);

//    auto **Verr = new double * [nx+1];
    double Verr[nx+6][nx+6] = {0.};

    //obliczamy blad
    for(int i=1; i<nx; i++){
        for(int j=1; j<ny; j++){
            Verr[i][j] = (Vn[i+1][j] - 2.*Vn[i][j] + Vn[i-1][j]) / (delta*delta) + (Vn[i][j+1] - 2.*Vn[i][j] + Vn[i][j-1])/(delta*delta) + P[i][j]/eps;
        }
    }

    //zapisujemy rzowiazania V i blad do pliku
    for(int i=1; i<nx; i++){
        for(int j=1; j<ny; j++){
            fprintf(fp_V,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, Vn[i][j]);
            fprintf(fp_Err,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, Verr[i][j]);
        }
    }
    std::cout<<"... Done\n";
}

void lokalna(double wl, std::string filename_S){
    std::cout<<"Relaksacja lokalna dla wl = " + std::to_string(wl) + " ...\n";
    FILE *fp_S = fopen(filename_S.c_str(),"w");
    double V [nx+5][ny+5] = {0.};
    double P [nx+5][ny+5] = {0.};
    //tablica gestosci i wstepne rozwiazania Vs Vn
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            V[i][0] = V1;
            P[i][j] = gestosc1(i*delta, j*delta) + gestosc2(i*delta, j*delta);
        }
    }

    int count = 0;
    double S[100000] = {0.};
    S[0] = 0.;

    do{
        count++;

        for(int i=1; i<nx; i++){
            for(int j=1; j<ny; j++){
                V[i][j] = (1-wl) * V[i][j] + (wl*0.25)*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + (delta*delta)/eps * P[i][j]);
            }
        }

        for(int j=1; j<ny; j++){
            V[0][j] = V[1][j];
            V[nx][j] = V[nx-1][j];
        }

        //warunek stopu
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                S[count] += (delta*delta) * (0.5 * std::pow((V[i+1][j] - V[i][j])/delta, 2) + 0.5*std::pow((V[i][j+1]-V[i][j])/delta, 2 ) - (P[i][j]*V[i][j]) );
            }
        }

        fprintf(fp_S,"%d\t%12.6g\n",count-1,S[count]);

    }while(std::abs((S[count] - S[count-1]) / S[count-1]) > TOL);
    std::cout<<"... Done\n";
}