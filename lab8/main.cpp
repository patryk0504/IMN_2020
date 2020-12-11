#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <omp.h>

const int nx = 400;
const int ny = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 10*delta;
const double xa = 0.45;
const double ya = 0.45;
const int IT_MAX = 11000;

void fillPolePredkosci(std::vector<std::vector<double>> &Vx, std::vector<std::vector<double>> &Vy, std::vector<std::vector<double>> &psi){
    for(int i=1; i<nx; i++){
        for(int j=1; j<ny; j++){
            Vx[i][j] = (psi[i][j+1] - psi[i][j-1])/(2.*delta);
            Vy[i][j] = - (psi[i+1][j] - psi[i-1][j])/(2.*delta);
        }
    }
    //na zatawce ustalamy
    for(int i = i_1; i <= i_2; i++){
        for(int j = 0; j<= j_1; j++){
            Vx[i][j] = 0.;
            Vy[i][j] = 0.;
        }
    }
    //na dolnym i gornym brzegu
    for(int i = 1; i<=nx-1; i++){
        Vx[i][0] = 0.;
        Vy[i][ny] = 0.;
    }
    //na lewym i prawym brzegu przepisujemy z sasiednich wezlow
    for(int j = 0; j<=ny; j++){
        Vx[0][j] = Vx[1][j];
    }
    for(int j=0; j<=ny; j++) {
        Vx[nx][j] = Vx[nx - 1][j];
    }
}

double modulV(double Vx, double Vy){
    return std::sqrt(std::pow(Vx,2) + std::pow(Vy,2));
}

double maxV(std::vector<std::vector<double>> &Vx, std::vector<std::vector<double>> &Vy){
    double max = modulV(Vx[0][0], Vy[0][0]);
    int iMax = 0;
    int jMax = 0;
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            if(modulV(Vx[i][j], Vy[i][j]) > max){
                max = modulV(Vx[i][j], Vy[i][j]);
                iMax = i;
                jMax = j;
            }
        }
    }

    return modulV(Vx[iMax][jMax], Vy[iMax][jMax]);
}

bool zastawka(int i, int j){
    if(i >=i_1 && i <= i_2 && j >= 0 && j <= j_1){
        return true;
    }else
        return false;
}

void algorytm(double D, double dt, std::vector<std::vector<double>> &Vx,std::vector<std::vector<double>> &Vy,
              std::vector<std::vector<double>> &U0, std::vector<std::vector<double>> &U1, FILE *result, std::string filename){
    for(int it = 0; it<= IT_MAX; it++){
        //start iteracji picarda
        for(int i=0; i<=nx; i++){     // U1 = U0
            for(int j=0; j<=ny; j++){
                U1[i][j] = U0[i][j];
            }
        }
        #pragma omp parallel for num_threads(6) collapse(3)
        for(int k=1; k<=20; k++){
            for(int i=0; i<=nx; i++){
                for(int j=1; j<=ny-1; j++){
                    if(zastawka(i,j)) {
                        continue;
                        //wzor 9 + periodyczne WB
                    }else if(i == 0){//warunek brzegowy (i = i-1) => (i = nx)
                        U1[i][j] = (1. / (1 + (2.*D*dt)/(delta*delta))) * ( U0[i][j] - ((dt/2.) * Vx[i][j]) * (  (U0[i+1][j] - U0[nx][j])/(2*delta) + (U1[i+1][j] - U1[nx][j])/(2*delta) )
                                                                            - ((dt/2)*Vy[i][j]) * (((U0[i][j+1] - U0[i][j-1])/(2*delta)) + ((U1[i][j+1] - U1[i][j-1])/(2*delta)))
                                                                            + ((dt/2)*D) * ((U0[i+1][j] + U0[nx][j] + U0[i][j+1] + U0[i][j-1] - 4.*U0[i][j])/(delta*delta) + (U1[i+1][j] + U1[nx][j] + U1[i][j+1] + U1[i][j-1])/(delta*delta))
                        );
                    }else if(i == nx){//wb (i = i+1) => (i = 0)
                        U1[i][j] = (1. / (1 + (2*D*dt)/(delta*delta))) * ( U0[i][j] - ((dt/2) * Vx[i][j]) * (  (U0[0][j] - U0[i-1][j])/(2*delta) + (U1[0][j] - U1[i-1][j])/(2*delta) )
                                                                           - ((dt/2)*Vy[i][j]) * (((U0[i][j+1] - U0[i][j-1])/(2*delta)) + ((U1[i][j+1] - U1[i][j-1])/(2*delta)))
                                                                           + ((dt/2)*D) * ((U0[0][j] + U0[i-1][j] + U0[i][j+1] + U0[i][j-1] - 4*U0[i][j])/(delta*delta) + (U1[0][j] + U1[i-1][j] + U1[i][j+1] + U1[i][j-1])/(delta*delta))
                        );
                    }
                    else{
                        U1[i][j] = (1. / (1 + (2*D*dt)/(delta*delta))) * ( U0[i][j] - ((dt/2) * Vx[i][j]) * (  (U0[i+1][j] - U0[i-1][j])/(2*delta) + (U1[i+1][j] - U1[i-1][j])/(2*delta) )
                                                                           - ((dt/2)*Vy[i][j]) * (((U0[i][j+1] - U0[i][j-1])/(2*delta)) + ((U1[i][j+1] - U1[i][j-1])/(2*delta)))
                                                                           + ((dt/2)*D) * ((U0[i+1][j] + U0[i-1][j] + U0[i][j+1] + U0[i][j-1] - 4*U0[i][j])/(delta*delta) + (U1[i+1][j] + U1[i-1][j] + U1[i][j+1] + U1[i][j-1])/(delta*delta))
                        );
                    }
                }
            }
        }

        for(int k = 0; k<50; k++){
            int T = k*IT_MAX/50;
            if(T == it){
                std::string filename2 = filename + std::to_string(k) + ".dat";
                FILE *map1 = fopen(filename2.c_str(),"w");
                for(int i=0; i<=nx; i++){
                    for(int j=0; j<=ny; j++){
                        fprintf(map1,"%d\t%d\t%12.6g\n",i,j,U1[i][j]);
                    }
                    fprintf(map1,"\n");
                }
                fclose(map1);
            }
        }

        //zachowujemy rozwiazanie do nast wywolania
        for( int i=0; i<=nx; i++ ){
            for( int j=0; j<=ny; j++ ) {
                U0[i][j] = U1[i][j];
            }
        }

        //wyznaczamy i zapisujemy do pliku c oraz xsr
        double c=0.;
        double xsr=0.;
        for( int i=0; i<=nx; i++ ){
            for( int j=0; j<=ny; j++ ){
                c += U0[i][j] * delta * delta;
                xsr += i*delta * U0[i][j] * delta * delta;
            }
        }

        fprintf(result,"%12.6g\t%12.6g\t%12.6g\n",it*dt, c,xsr);
        printf("iteracja: %d c: %f xsr: %f\n",it,c,xsr);
    }
}

int main() {

    FILE *result = fopen("res.dat","w");
    FILE *resultVx = fopen("resVx.dat","w");
    FILE *resultVy = fopen("resVy.dat","w");

    FILE *result2 = fopen("res2.dat","w");

    std::vector<std::vector<double>> psi;
    psi.resize(nx+1, std::vector<double>(ny + 1, 0.0));

    std::vector<std::vector<double>> Vx;
    Vx.resize(nx+1, std::vector<double>(ny + 1, 0.0));

    std::vector<std::vector<double>> Vy;
    Vy.resize(nx+1, std::vector<double>(ny + 1, 0.0));

    std::vector<std::vector<double>> U0;
    U0.resize(nx+1, std::vector<double>(ny + 1, 0.0));

    std::vector<std::vector<double>> U1;
    U1.resize(nx+1, std::vector<double>(ny + 1, 0.0));

    std::ifstream infile("psi.dat");
    int i;
    int j;
    double psi_v;
    while(infile >> i >> j >> psi_v){
        psi[i][j] = psi_v;
    }
    i=j=0;

    fillPolePredkosci(Vx,Vy,psi);
    //zapis Vx i Vy
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            fprintf(resultVx,"%d\t%d\t%12.6g\n",i,j,Vx[i][j]);
            fprintf(resultVy,"%d\t%d\t%12.6g\n",i,j,Vy[i][j]);
        }
        fprintf(resultVx,"\n");
        fprintf(resultVy,"\n");
    }
    double Vmax = maxV(Vx,Vy);
    double dt = delta/(4*Vmax);
    /////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    double D = 0;
    std::cout<<"Vmax: "<<Vmax<<" dt: "<<dt<<std::endl;
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            U0[i][j] = 1./(2.*M_PI*sigma*sigma) * 1.* std::exp(-1. * (std::pow(i*delta - xa,2) + std::pow(j*delta-ya,2))/(2.*sigma*sigma));
        }
    }
    std::string filename = "mapa_";
    algorytm(D,dt,Vx,Vy,U0,U1,result,filename);

    //////////////////////////////////////////////////////////////
    D = 0.1;
    for(int i=0; i<=nx; i++){//wracamy do warunkow poczatkowych
        for(int j=0; j<=ny; j++){
            U0[i][j] = 1./(2.*M_PI*sigma*sigma) * 1.* std::exp(-1. * (std::pow(i*delta - xa,2) + std::pow(j*delta-ya,2))/(2.*sigma*sigma));
        }
    }
    filename = "mapa2_";
    algorytm(D,dt,Vx,Vy,U0,U1,result2,filename);
    //////////////////////////////////////////////////////////////

    fclose(result);
    fclose(resultVx);
    fclose(resultVy);
    fclose(result2);

    return 0;
}
