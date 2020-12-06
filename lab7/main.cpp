#include <iostream>
#include <stdio.h>
#include <cmath>
#include "vector"

//////////////////////////////////////// parametry zadania
const double delta = 0.01;
const int nx = 200;
const int ny = 90;
const int IT_MAX = 20000;
const int i_1 = 50;
const int j_1 = 55;
const int u = 1.;
const double p = 1.;
////////////////////////////////////////

double Qwy(double Qwe, double yny, double yj1){//zwraca Qwyjsciowe na podstawie wejsciowego
    return Qwe * (std::pow(yny,3) - std::pow(yj1,3) - 3*yj1*std::pow(yny,2) + 3*std::pow(yj1,2)*yny)/std::pow(yny,3);
}
bool checkBrzeg(int i, int j){//true - indeksy nie odnosza sie do brzegu
    //sprawdzam krawedzie
    return (i!=0 && j !=0 && i!=nx && j!=ny && !(i<=i_1 && j==j_1) && !(i==i_1 && j<=j_1) && !(i<=i_1 && j < j_1));//true - nie na krawedzi
}
double wzor_8(std::vector<std::vector<double>> &fi,std::vector<std::vector<double>> &C, int i, int j){
    return 0.25 * (fi[i+1][j] + fi[i-1][j] + fi[i][j+1] + fi[i][j-1] - delta*delta*C[i][j]);
}
double wzor_9(std::vector<std::vector<double>> &fi, std::vector<std::vector<double>> &C, int i, int j, int omega){
    return 0.25 * (C[i+1][j] + C[i-1][j] + C[i][j+1] + C[i][j-1]) - omega*p/(16.0*u)*(((fi[i][j+1] - fi[i][j-1]) * (C[i+1][j] - C[i-1][j]) - (fi[i+1][j] - fi[i-1][j]) * (C[i][j+1] - C[i][j-1])));
}
//funkcje aktualizujace macierze psi i ksi
auto WB_1(double Qwe, std::vector<std::vector<double>> &fi);
auto WB_2(double Qwe, std::vector<std::vector<double>> &fi,std::vector<std::vector<double>> &C);
//glowna funkcja przeprowadzajaca relaksacje
void relaksacja(double Qwe, const std::string& filename_kontury, const std::string& filename_predkosci);



int main() {
//zad1
    double Qwe = -1000;
    relaksacja(Qwe,"fi_c_minus1000.dat","u_v_minus1000.dat");
//zad2
    Qwe = -4000;
    relaksacja(Qwe, "fi_c_minus4000.dat","u_v_minus4000.dat");
//zad3
    Qwe = 4000;
    relaksacja(Qwe,"fi_c_plus4000.dat","u_v_plus4000.dat");
    return 0;
}

auto WB_1(double Qwe, std::vector<std::vector<double>> &fi){
    double yj1 = delta*j_1;
    double yny = delta*ny;
    //brzeg A - wejscie
    for(int j = j_1; j<ny+1; j++){
        double y = delta*j;
        fi[0][j] = Qwe/(2*u) * (std::pow(y,3)/3 - std::pow(y,2)/2 * (yj1 + yny)+y*yj1*yny);
    }
    //brzeg B
    for(int i = 1; i< nx; i++){
        fi[i][ny] = fi[0][ny];
    }
    //brzeg C
    for(int j=0; j<ny+1; j++){
        double y = j*delta;
        fi[nx][j] = Qwy(Qwe,yny,yj1)/(2*u) * (std::pow(y,3)/3 - std::pow(y,2)/2 * yny) + (Qwe*std::pow(yj1,2)*(-1.*yj1+3*yny))/(12*u);
    }
    //brzeg D
    for(int i = i_1; i<nx; i++){
        fi[i][0] = fi[0][j_1];
    }
    //brzeg E
    for(int j = 1; j<j_1+1; j++){
        fi[i_1][j] = fi[0][j_1];
    }
    //brzeg F
    for(int i = 1; i<i_1+1; i++){
        fi[i][j_1] = fi[0][j_1];
    }
}

auto WB_2(double Qwe, std::vector<std::vector<double>> &fi,std::vector<std::vector<double>> &C){
    double yj1 = delta*j_1;
    double yny = delta*ny;
    //brzeg A - wejscie
    for(int j = j_1; j<ny+1; j++){
        double y = delta*j;
        C[0][j] = Qwe/(2*u) * (2*y - yj1 - yny);
    }
    //brzeg C - wyjscie
    for(int j = 0; j< ny+1; j++){
        double y = delta * j;
        C[nx][j] = Qwy(Qwe, yny, yj1) / (2*u) * (2*y - yny);
    }
    //brzeg B
    for(int i = 1; i<nx; i++){
        C[i][ny] = 2/(delta*delta) * (fi[i][ny-1] - fi[i][ny]);
    }
    //brzeg D
    for(int i=i_1+1; i<nx;i++){
        C[i][0] = 2/(delta*delta) * (fi[i][1] - fi[i][0]);
    }
    //brzeg E
    for(int j=1; j<j_1+1; j++){
        C[i_1][j] = 2/(delta*delta) * (fi[i_1+1][j]-fi[i_1][j]);
    }
    //brzeg F
    for(int i = 1; i<i_1+1; i++){
        C[i][j_1] = 2/(delta*delta) * (fi[i][j_1+1] - fi[i][j_1]);
    }
    //wierzcholek E/F
    C[i_1][j_1] = 0.5 * (C[i_1-1][j_1] + C[i_1][j_1-1]);
}

void relaksacja(double Qwe, const std::string& filename_kontury, const std::string& filename_predkosci){
    FILE * kont = fopen(filename_kontury.c_str(),"w");
    FILE * pred = fopen(filename_predkosci.c_str(),"w");

    std::vector<std::vector<double>> fi;
    fi.resize(nx+1, std::vector<double>(ny + 1, 0));
    std::vector<std::vector<double>> C;
    C.resize(nx+1, std::vector<double>(ny + 1, 0));
    WB_1(Qwe,fi);
    int omega;

    for(int it = 1; it<IT_MAX; it++){
        if(it < 2000)
            omega = 0;
        else
            omega = 1;
        for(int i = 1; i<nx; i++){
            for(int j=1; j<ny; j++){
                if(checkBrzeg(i,j)){
                    fi[i][j] = wzor_8(fi,C,i,j);
                    C[i][j] = wzor_9(fi,C,i,j,omega);
                }
            }
        }
        WB_2(Qwe,fi,C);
        double blad = 0;
        for(int i=1; i<nx; i++){
            blad += (fi[i+1][j_1 + 2] + fi[i-1][j_1 + 2] + fi[i][j_1 + 2 + 1] + fi[i][j_1+2 -1] - 4*fi[i][j_1+2] - delta*delta*C[i][j_1+2]);
        }
        std::cout<<"iteracja: "<<it<<" blad: "<<blad<<std::endl;
    }
    //fi i c
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            fprintf(kont,"%12.6g\t%12.6g\t%12.6g\t%12.6g\n",i*delta,j*delta,fi[i][j],C[i][j]);
        }
        fprintf(kont,"\n");
    }
    std::vector<std::vector<double>> u_pred_vec;
    u_pred_vec.resize(nx+1, std::vector<double>(ny + 1, 0));
    std::vector<std::vector<double>> v_pred_vec;
    v_pred_vec.resize(nx+1, std::vector<double>(ny + 1, 0));

    double u_pred;
    double v_pred;
//u i v
    for(int i=1; i<nx; i++){
        for(int j=1; j<ny; j++){
            if(checkBrzeg(i,j)) {
                u_pred = (fi[i][j + 1] - fi[i][j - 1]) / (2 * delta);
                v_pred = -(fi[i + 1][j] - fi[i - 1][j]) / (2 * delta);

                u_pred_vec[i][j] = u_pred;
                v_pred_vec[i][j] = v_pred;
            }
        }
    }
    for(int i = 1; i< nx; i++){
        for(int j = 1; j<ny; j++){
            fprintf(pred,"%12.6g %12.6g %f %f\n",i*delta,j*delta,u_pred_vec[i][j],v_pred_vec[i][j]);
        }
        fprintf(pred,"\n");
    }

    fclose(kont);
    fclose(pred);
}