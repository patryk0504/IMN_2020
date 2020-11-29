#include <iostream>
#include "mgmres.h"
//#include "mgmres.c"
#include <stdio.h>
#include <cmath>
#include <array>

//niektore parametry wejsciowe
const double delta = 0.1;
int V1 = 10;
int V2 = -10;
int V3 = 10;
int V4 = -10;
//funkcje zwracajace odpowiednie parametry dla zmiennych parametrow wejsciowych
double xmax(int nx){return delta*nx;}
double ymax(int ny){return delta*ny;}
int N(int nx, int ny){return (nx+1)*(ny+1);}
int j(int nx, int l){return (l/(nx+1));}
int i(int nx, int l){return (l - j(nx,l)*(nx+1));}
double element(int nx, int l, int eps1, int eps2){return i(nx,l) <= nx/2 ? eps1 : eps2;}

double gesetosc1(double x, double y, int nx, int ny, double xmax, double ymax){
    double sigma = xmax/10.;
    return std::exp(-1. * std::pow(x-0.25*xmax,2)/(sigma*sigma) - std::pow(y - 0.5*ymax,2)/(sigma*sigma));
}
double gestosc2(double x, double y, int nx, int ny, double xmax, double ymax){
    double sigma = xmax/10.;
    return -1. * std::exp(-1. * std::pow(x-0.75*xmax,2)/(sigma*sigma) - std::pow(y - 0.5*ymax,2)/(sigma*sigma));
}
//funkcja odpowiadajaca za wypelnienie wektorow i macierzy
int wypelnienieDirichl(int *ja, int *ia, double *a, double *b,double xmax,double ymax,int nx, int ny, int eps1, int eps2,  FILE * macierz, FILE * wektor, bool flag = false);
//funkcja uruchamijaca procedure z biblioteki mgmres z odpowiednimi parametrami
void algorytm(int nx, int ny, bool flag = true, std::string filename = "none", int eps1 = 1, int eps2 = 1, double xmax = 0, double ymax=0);

int main() {

    int nx = 4;
    int ny = 4;
    int eps1 = 1;
    int eps2 = 1;
    //stale wartosci eps1 i eps2 ->
    //wyznaczenie wektora b i macierzy A
    algorytm(nx,ny);
    //mapa potencialu a
    nx = 50;
    ny = 50;
    algorytm(nx,ny,false,"map50_50.txt");
    //mapa potencialu b
    nx = 100;
    ny = 100;
    algorytm(nx,ny,false, "map100_100.txt");
    //mapa potencialu c
    nx = 200;
    ny = 200;
    algorytm(nx,ny,false, "map200_200.txt");
    //zad6
    nx = 100;
    ny = 100;
    V1 = V2 = V3 = V4 = 0;
    //bedziemy zmieniac eps1 i eps2 ->
    //a
    eps1 = 1;
    eps2 = 1;
    algorytm(nx,ny,false, "map_e_1_1.txt",eps1, eps2,xmax(nx),ymax(ny));
    //b
    eps1 = 1;
    eps2 = 2;
    algorytm(nx,ny, false,"map_e_1_2.txt",eps1, eps2, xmax(nx), ymax(ny));
    //c
    eps1 = 1;
    eps2 = 10;
    algorytm(nx,ny, false,"map_e_1_10.txt",eps1, eps2, xmax(nx), ymax(ny));

    return 0;
}

int wypelnienieDirichl(int *ja, int *ia, double *a, double *b,double xmax, double ymax, int nx, int ny, int eps1, int eps2, FILE * macierz, FILE * wektor, bool flag){
    //numeruje niezerowe elem A
    int k = -1;
    //liczba niezerowych el (1 el ma indeks 0)
    for(int l = 0; l < N(nx, ny); ++l) {
        int brzeg = 0;   // wskaźnik położenia : 0 - środek obszaru ; 1 - brzeg
        double vb = 0;   // potencjal na brzegu
        if(i(nx, l) == 0){// lewy brzeg
            brzeg = 1;
            vb = V1;
        }
        if(j(nx, l) == ny){// górny brzeg
            brzeg = 1;
            vb = V2;
        }
        if(i(nx, l) == nx){// prawy brzeg
            brzeg = 1;
            vb = V3;
        }
        if(j(nx, l) == 0){// dolny brzeg
            brzeg = 1;
            vb = V4;
        }
        // wypełniamy od razu wektor wyrazów wolnych
        b[l] = (-1)*(gesetosc1(delta*i(nx,l), delta*j(nx,l),nx,ny,xmax,ymax) + gestosc2(delta*i(nx,l), delta*j(nx,l),nx,ny,xmax,ymax)); //sigma

        if(brzeg == 1)
            b[l] = vb; // wymuszamy wartość potencjału na brzegu

        ia[l] = -1;// wskaźnik dla pierwszego el . w wierszu
        // lewa skrajna przekatna
        if(l - nx - 1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0)
                ia[l] = k;

            a[k] = element(nx,l,eps1,eps2)/(delta*delta);
            if(flag)    fprintf(macierz,"%d %0.f \n", k, a[k]);
            ja[k] = l - nx - 1;
        }
        // poddiagonala
        if(l-1 >= 0 && brzeg == 0) {
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = element(nx,l,eps1,eps2)/(delta*delta);
            if(flag)    fprintf(macierz,"%d %0.f \n", k, a[k]);
            ja[k] = l-1;
        }
        //diagonala
        k++;
        if(ia[l] < 0)
            ia[l] = k;

        if(brzeg == 0) {
            a[k] = -(2 * element(nx, l, eps1, eps2) + element(nx, l + 1, eps1, eps2) +
                     element(nx, l + nx + 1, eps1, eps2)) / (delta * delta);
            if(flag)    fprintf(macierz,"%d %0.f \n", k, a[k]);
        }else {
            a[k] = 1;
            if(flag)    fprintf(macierz,"%d %0.f \n", k, a[k]);
        }
        ja[k] = l;
        //naddiagonala
        if(l < N(nx, ny) && brzeg == 0){
            k++;
            a[k] = element(nx,l+1,eps1,eps2)/(delta*delta);
            if(flag)    fprintf(macierz,"%d %0.f \n", k, a[k]);
            ja[k] = l + 1;
        }
        //prawa skrajna przekatna
        if(l < N(nx, ny)-nx-1 && brzeg == 0){
            k++;
            a[k] = element(nx,l+nx+1,eps1,eps2)/(delta*delta);
            if(flag)    fprintf(macierz,"%d %0.f \n", k, a[k]);
            ja[k] = l + nx + 1;
        }
        if(flag){
            if(l%5 == 0 && l != 0)
                fprintf(wektor, "\n");
            fprintf(wektor,"%d %d %d %f \n", l, i(nx,l), j(nx,l), b[l]);
        }
    }
    int nz_num = k+1;
    ia[N(nx, ny)] = nz_num;
    if(flag) {
        for (int z = nz_num; z < 5 * N(nx, ny); z++)
            fprintf(macierz, "%d %d \n", z, 0);
    }
    return nz_num;
}

void algorytm(int nx, int ny, bool flag, std::string filename,int eps1, int eps2, double xmax, double ymax){
    if(flag){
        FILE *resMacierz = fopen("macierz.txt","w");
        FILE *resWektor = fopen("wektor.txt","w");
        double a[5*N(nx,ny)];
        int ja[5*N(nx,ny)];
        int ia[N(nx,ny)+1];
        double b [N(nx,ny)];
        int nz_num = wypelnienieDirichl(ja,ia,a,b,xmax,ymax,nx,ny,eps1,eps2,resMacierz,resWektor, true);
        fclose(resWektor);
        fclose(resMacierz);
    }else{
        FILE *resMapa = fopen(filename.c_str(), "w");
        auto *a = new double [5*N(nx,ny)];
        int *ja = new int [5*N(nx,ny)];
        int *ia = new int[N(nx,ny)+1];
        auto *b = new double [N(nx,ny)];
        auto *V = new double [N(nx,ny)];
        ///parametry do procedury pmgmres_ilu_cr
        int nz_num = wypelnienieDirichl(ja,ia,a,b,xmax,ymax, nx,ny,eps1,eps2,nullptr, nullptr);
        int itr_max = 500;
        int mr = 500;
        double tol_abs = pow(10,-8);
        double tol_rel = pow(10,-8);
        /////////////////////////////////////////////////////
        pmgmres_ilu_cr(N(nx,ny),nz_num,ia,ja,a,V,b,itr_max, mr,tol_abs,tol_rel);
        double max = 0;
        for(int k=0; k<N(nx,ny); ++k){
            if(delta * i(nx,k) < max){// aby uzyc fukcji splot gnuplot checmy miec puste wiersze pomiedzy odpowiednimi rekordami
                fprintf(resMapa,"\n");
            }
            fprintf(resMapa,"%f %f %f \n", delta*i(nx,k),delta*j(nx,k), V[k]);
            max = delta * i(nx,k);
        }
        fclose(resMapa);
        free(a);
        free(ja);
        free(ia);
        free(b);
        free(V);
    }
}
