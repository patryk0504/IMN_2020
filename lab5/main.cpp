#include <iostream>
#include <cmath>
#include <array>
#include <stdio.h>

//przyjete wartosci parametrow
double delta = 0.2;
const int nx = 128;
const int ny = 128;
const int nx1 = 129;
const int ny1 = 129;
double xmax = delta*nx;
double ymax = delta*ny;
double TOL = std::pow(10,-8);
int iteracje = 0;
/////////////////////////////////////////////////


//wywolanie dla kazdego k = 16,8,4,2,1
void relWielosiatkowa(int k, std::array<std::array<double, nx1>,ny1> &arr, std::string filename_s, std::string filename_v, bool flag);

//funkcja przygotowujaca tablice V
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

    return 0;
}

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
                arr[i+k][j+k/2] = i == nx - k ? arr[i+k][j+k/2] : 0.5 * (arr[i+k][j] + arr[i+k][j+k]);
                arr[i+k/2][j+k] = j == ny - k ? arr[i+k/2][j+k] : 0.5 * (arr[i][j+k] + arr[i+k][j+k]);
                arr[i+k/2][j] = j == 0 ? arr[i+k/2][j] : 0.5 * (arr[i][j] + arr[i+k][j]);
                arr[i][j+k/2] = i == 0 ? arr[i][j+k/2] : 0.5 * (arr[i][j] + arr[i][j+k]);
            }
        }
    }
    fclose(sFile);
    fclose(vFile);
}
