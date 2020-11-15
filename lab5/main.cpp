#include <iostream>
#include <cmath>
#include <array>
#include <stdio.h>

double delta = 0.2;
const int nx = 128;
const int ny = 128;
const int nx1 = 129;
const int ny1 = 129;
double xmax = delta*nx;
double ymax = delta*ny;
double TOL = std::pow(10,-8);
int iteracje = 0;

//wersja wywolanie dla kazdego k = 16,8,4,2,1
void relWielosiatkowa(int k, std::array<std::array<double, nx1>,ny1> &arr, std::string filename_s, std::string filename_v, bool flag){
    FILE *sFile = fopen(filename_s.c_str(), "w");
    FILE *vFile = fopen(filename_v.c_str(), "w");
    double s[100000] ={0.};
    do{
        iteracje++;
        for(int i = k; i<=nx-k; i+=k){
            for(int j = k; j<=ny-k; j+=k){
                arr[i][j] = 0.25 * (arr[i+k][j] + arr[i-k][j] + arr[i][j+k] + arr[i][j-k]);
            }
        }

        for(int i=0; i<=nx-k; i+=k){
            for(int j=0; j<=ny-k; j+=k){
                s[iteracje] += std::pow(k*delta,2)/2. * (std::pow((arr[i+k][j] - arr[i][j]) / (2*k*delta) + (arr[i+k][j+k] - arr[i][j+k])/(2*k*delta),2) +
                        std::pow((arr[i][j+k] - arr[i][j])/(2*k*delta) + (arr[i+k][j+k] - arr[i+k][j])/(2*k*delta),2));
            }
        }
        fprintf(sFile,"%d\t%12.6g\n",iteracje-1,s[iteracje]);
    }while(std::abs((s[iteracje]-s[iteracje-1])/s[iteracje-1]) > TOL);

    std::cout<<"k= "<<k<<" iteracje: " << iteracje<<std::endl;

    for(int i = 0; i<= nx; i+=k){
        for(int j = 0; j<= ny; j+=k){
            fprintf(vFile,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, arr[i][j]);
        }
    }

    if(flag){
        for(int i=0; i<=nx-k; i+=k){
            for(int j=0; j<=ny-k; j+=k){
                arr[i+k/2][j+k/2] = 0.25 * (arr[i][j] + arr[i+k][j] + arr[i][j+k] + arr[i+k][j+k]);
                if(i!=nx-k)
                    arr[i+k][j+k/2] = 0.5 * (arr[i+k][j] + arr[i+k][j+k]);
                if(j!=ny-k)
                    arr[i+k/2][j+k] = 0.5 * (arr[i][j+k] + arr[i+k][j+k]);
                if(j!=0)
                    arr[i+k/2][j] = 0.5 * (arr[i][j] + arr[i+k][j]);
                if(i!=0)
                    arr[i][j+k/2] = 0.5 * (arr[i][j] + arr[i][j+k]);
            }
        }
    }
}



/* // wersja wszystko w jednej metodzie
void relWielosiatkowa2(std::array<std::array<double, nx1>, ny1> &arr){
    std::string fileS[]={"S_16.txt","S_8.txt","S_4.txt","S_2.txt","S_1.txt"};
    std::string fileV[]={"V_16.txt","V_8.txt","V_4.txt","V_2.txt","V_1.txt"};
    FILE *sFile16 = fopen(fileS[0].c_str(), "w");
    FILE *vFile16 = fopen(fileV[0].c_str(), "w");
    FILE *sFile8 = fopen(fileS[1].c_str(), "w");
    FILE *vFile8 = fopen(fileV[1].c_str(), "w");
    FILE *sFile4 = fopen(fileS[2].c_str(), "w");
    FILE *vFile4 = fopen(fileV[2].c_str(), "w");
    FILE *sFile2 = fopen(fileS[3].c_str(), "w");
    FILE *vFile2 = fopen(fileV[3].c_str(), "w");
    FILE *sFile1 = fopen(fileS[4].c_str(), "w");
    FILE *vFile1 = fopen(fileV[4].c_str(), "w");

    int iteracje = 0;

    for(int k = 16; k>0; k/=2){
        double s[100000] ={0.};
//        int iteracje = 0;
        do{
            iteracje++;
            for(int i = k; i<=nx-k; i+=k){
                for(int j = k; j<=ny-k; j+=k){
                    arr[i][j] = 0.25 * (arr[i+k][j] + arr[i-k][j] + arr[i][j+k] + arr[i][j-k]);
                }
            }

            for(int i=0; i<=nx-k; i+=k){
                for(int j=0; j<=ny-k; j+=k){
                    s[iteracje] += std::pow(k*delta,2)/2. * (std::pow((arr[i+k][j] - arr[i][j]) / (2*k*delta) + (arr[i+k][j+k] - arr[i][j+k])/(2*k*delta),2) +
                                                             std::pow((arr[i][j+k] - arr[i][j])/(2*k*delta) + (arr[i+k][j+k] - arr[i+k][j])/(2*k*delta),2));
                }
            }

//        std::cout<<iteracje-1<<" "<<s[iteracje]<<std::endl;
        if(k==16) {
            fprintf(sFile16, "%d\t%12.6g\n", iteracje - 1, s[iteracje]);
        }else if(k==8){
            fprintf(sFile8, "%d\t%12.6g\n", iteracje - 1, s[iteracje]);

        }else if(k==4){
            fprintf(sFile4, "%d\t%12.6g\n", iteracje - 1, s[iteracje]);

        }else if(k==2){
            fprintf(sFile2, "%d\t%12.6g\n", iteracje - 1, s[iteracje]);

        }else if(k==1){
            fprintf(sFile1, "%d\t%12.6g\n", iteracje - 1, s[iteracje]);

        }

        }while(std::abs((s[iteracje]-s[iteracje-1])/s[iteracje-1]) > TOL);

        std::cout<<"k= "<<k<<" iteracje: " << iteracje<<std::endl;

        for(int i = 0; i<= nx; i+=k){
            for(int j = 0; j<= ny; j+=k){
                if(k==16)
                    fprintf(vFile16,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, arr[i][j]);
                else if(k==8)
                    fprintf(vFile8,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, arr[i][j]);
                else if(k==4)
                    fprintf(vFile4,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, arr[i][j]);
                else if(k==2)
                    fprintf(vFile2,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, arr[i][j]);
                else if(k==1)
                    fprintf(vFile1,"%12.6g\t%12.6g\t%12.6g\n",i*delta, j*delta, arr[i][j]);

            }
        }

        for(int i=0; i<=nx-k; i+=k){
            for(int j=0; j<=ny-k; j+=k){
                arr[i+k/2][j+k/2] = 0.25 * (arr[i][j] + arr[i+k][j] + arr[i][j+k] + arr[i+k][j+k]);
                if(i!=nx-k)
                    arr[i+k][j+k/2] = 0.5 * (arr[i+k][j] + arr[i+k][j+k]);
                if(j!=ny-k)
                    arr[i+k/2][j+k] = 0.5 * (arr[i][j+k] + arr[i+k][j+k]);
                if(j!=0)
                    arr[i+k/2][j] = 0.5 * (arr[i][j] + arr[i+k][j]);
                if(i!=0)
                    arr[i][j+k/2] = 0.5 * (arr[i][j] + arr[i][j+k]);
            }
        }

    }
}

 */

auto prepareV(){
    std::array<std::array<double, nx1>, ny1> arr = {{0.}};
    for(int x=0; x<=nx; x++){
        arr[x][ny] = -1 * std::sin(2*M_PI*delta*x/xmax);
        arr[x][0] =std::sin(2*M_PI*delta*x/xmax);
    }
    for(int y=0; y<=ny; y++){
        arr[0][y] =std::sin(M_PI*delta*y/ymax);
        arr[nx][y] =std::sin(M_PI*delta*y/ymax);
    }
    return arr;
}

int main() {
    auto arr = prepareV();

    relWielosiatkowa(16,arr,"S_16.txt","V_16.txt", true);
    relWielosiatkowa(8,arr,"S_8.txt", "V_8.txt",true);
    relWielosiatkowa(4,arr,"S_4.txt","V_4.txt",true);
    relWielosiatkowa(2,arr,"S_2.txt", "V_2.txt",true);
    relWielosiatkowa(1,arr,"S_1.txt", "V_1.txt", true);

//    relWielosiatkowa2(arr);

    return 0;
}
