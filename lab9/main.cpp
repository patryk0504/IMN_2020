#include <iostream>
#include <stdio.h>
#include</usr/include/gsl/gsl_math.h>
#include</usr/include/gsl/gsl_linalg.h>
#include</usr/include/gsl/gsl_blas.h>

const int nx = 40;
const int ny = 40;
const int N = (nx + 1)*(ny + 1);
const int delta = 1;
const int dt = 1;
const int Ta = 40;
const int Tb = 0;
const int Tc = 30;
const int Td = 0;
double kb = 0.1;
double kd = 0.6;
const int IT_MAX = 2000;

int l(int i, int j){
    return i+j*(nx+1);
}

int main(){
    //tworzymy macierze A, B i wektor c
    gsl_matrix * A = gsl_matrix_calloc(N,N);
    gsl_matrix * B = gsl_matrix_calloc(N,N);
    gsl_vector *c = gsl_vector_calloc(N);

    //wypelniamy zgodnie z wzorami

    //niezerowe elem macierzowe z WB
    //wnetrze obszaru
    for(int i = 1; i<= nx - 1; i++){
        for(int j = 1; j <= ny-1; j++){
            double w1 = dt/(2.*delta*delta);
            gsl_matrix_set(A, l(i,j), l(i,j)-nx-1, w1);
            gsl_matrix_set(A,l(i,j),l(i,j)-1,w1);
            gsl_matrix_set(A, l(i,j), l(i,j)+1,w1);
            gsl_matrix_set(A, l(i,j), l(i,j)+nx+1,w1);
            w1 = -1.0 * (2*dt)/(delta*delta) - 1;
            gsl_matrix_set(A, l(i,j), l(i,j), w1);
            w1 = -(dt)/(2.0*delta*delta);
            gsl_matrix_set(B, l(i,j), l(i,j)-nx-1, w1);
            gsl_matrix_set(B, l(i,j), l(i,j) - 1, w1);
            gsl_matrix_set(B, l(i,j), l(i,j) + 1, w1);
            gsl_matrix_set(B, l(i,j), l(i,j)+nx+1, w1);
            w1 = (2.0*dt)/(delta*delta) - 1;
            gsl_matrix_set(B, l(i,j), l(i, j), w1);

        }
    }
    //WB Dirichleta
    for(int i = 0; i<=nx; i+=nx){
        for(int j = 0; j<=ny; j++){
            gsl_matrix_set(A, l(i,j), l(i, j), 1);
            gsl_matrix_set(B, l(i, j), l(i, j), 1);
            gsl_vector_set(c, l(i,j), 0);
        }
    }
    //WB von Neumanna na gornym brzegu
    for(int i=0; i<=nx-1; i++){
        int j = ny;
        double w1 = -1.0/(kb*delta);
        gsl_matrix_set(A, l(i,j), l(i, j) - nx - 1, w1);
        w1 = 1 + 1.0/(kb*delta);
        gsl_matrix_set(A, l(i,j), l(i,j), w1);
        gsl_vector_set(c, l(i,j), Tb);
        int lIndex = l(i,j);
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                gsl_matrix_set(B, lIndex, l(i,j), 0);
            }
        }
    }
    //WB von Neumanna na dolnym brzegu
    for(int i=1; i<=nx-1; i++){
        int j = 0;
        double w1 = 1 + 1.0/(kd*delta);
        gsl_matrix_set(A, l(i,j), l(i, j), w1);
        w1 = -1.0/(kd*delta);
        gsl_matrix_set(A, l(i,j), l(i,j) + nx + 1, w1);
        gsl_vector_set(c, l(i,j), Td);
        int lIndex = l(i,j);
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                gsl_matrix_set(B, lIndex, l(i,j), 0);
            }
        }
    }

    //tworzymy wektor startowy T i narzucamy war poczatkowe
    gsl_vector * T = gsl_vector_calloc(N);
    int i = 0;
    for(int j =0; j<=ny; j++){
        gsl_vector_set(T, l(i,j), Ta);
    }
    i = nx;
    for(int j =0; j<=ny; j++){
        gsl_vector_set(T, l(i,j), Tc);
    }
    for(int i=1; i<=nx-1; i++){
        for(int j=0; j<=ny; j++){
            gsl_vector_set(T, l(i,j),0);
        }
    }

    //znalezc rozklad LU macierzy A przy pomocy GSL
    gsl_permutation *p=gsl_permutation_calloc(N);
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);// p - wektor pamietajacy zamiane wierszy, signum - liczba permutacji


    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *Ttemp = gsl_vector_calloc(N);
    //implmentacja algorytmu CN
    for(int it=0; it<=IT_MAX; it++){
        gsl_blas_dgemv(CblasNoTrans, 1, B, T, 0, d);//d = A*x
        gsl_blas_daxpy(1, c, d); // d = d + c

        //zapisuje wektor T
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                gsl_vector_set(Ttemp,l(i,j),gsl_vector_get(T,l(i,j)));
            }
        }
        //rozw ukl rownani LU: A*T = d
        gsl_linalg_LU_solve(A,p,d,T);
        //zapisuje wektor T^n
        //zapisuje wektor T
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                gsl_vector_set(Ttemp, l(i,j), (gsl_vector_get(T,l(i,j)) - gsl_vector_get(Ttemp,l(i,j)))/dt );
            }
        }
        std::cout<<"teracja: "<<it<<std::endl;
        if(it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000){
            //mapy rozklady T w pomieszczeniu
            std::string filename = "mapaT_" + std::to_string(it) + ".dat";
            FILE *mapaT = fopen(filename.c_str(), "w");
            std::cout<<"zapisuje do pliku\n";
            for(int i=0; i<=nx; i++){
                for(int j=0; j<=ny; j++){
                    fprintf(mapaT,"%d\t%d\t%12.6g\n",i,j,gsl_vector_get(T,l(i,j)));
                }
                fprintf(mapaT,"\n");
            }
            fclose(mapaT);

            //mapy rozkladu
            filename = "mapaDT_" + std::to_string(it) + ".dat";
            FILE *mapaDT = fopen(filename.c_str(), "w");
            for(int i=0; i<=nx; i++){
                for(int j=0; j<=ny; j++){
                    fprintf(mapaDT,"%d\t%d\t%12.6g\n",i,j,gsl_vector_get(Ttemp,l(i,j)));
                }
                fprintf(mapaDT,"\n");
            }
            fclose(mapaDT);
        } 
    }


    return 0;
}