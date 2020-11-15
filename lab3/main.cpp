#include <iostream>
#include <stdio.h>
#include <cmath>
#include <array>
#include <functional>

double funF(double v, double t = 0, double x = 0){return v;}
double funG(double v, double x, double alfa, double t = 0){return alfa*(1. - std::pow(x,2)) * v - x;}

//metoda trapezów
auto metTrapezow(double xn, double vn, double dt, double alfa){//return std::array<double,2>
    double delta = std::pow(10,-10);
    double resultX = xn;
    double resultV = vn;
    double DX = 0;
    double DV = 0;
    do{
        double a11 = 1.;
        double a12 = -dt/2.;
        double a21 = (-dt/2) * (-2*alfa*resultX*resultV -1.);
        double a22 = 1. - (dt/2.) * alfa * (1 - std::pow(resultX,2));

        double F = resultX - xn - (dt/2.) * (funF(vn) + funF(resultV));
        double G = resultV - vn - (dt/2.) * (funG(vn,xn,alfa) + funG(resultV,resultX,alfa));

        DX = (-F * a22 - (-G)*a12)/(a11 * a22 - a12 * a21);
        DV = (a11 * (-G) - a21*(-F))/(a11*a22 - a12*a21);

        resultX += DX;
        resultV += DV;

    }while(DX >= delta || DV >= delta);

    return std::array<double,2>{resultX,resultV};
}
//metoda RK2
auto metRK2(double xn, double vn, double dt, double alfa){

    double k1x = funF(vn);
    double k1v = funG(vn,xn,alfa);

    double k2x = vn + dt*k1v;
    double k2v = funG(vn + dt*k1v,xn + dt*k1x,alfa);


    double resultX = xn + (dt/2.) * (k1x + k2x);
    double resultV = vn + (dt/2.) * (k1v + k2v);

    return std::array<double,2>{resultX,resultV};
}

void kontrola_kroku(double TOL, std::string filename, std::function<std::array<double,2>(double,double,double,double)> metoda);

int main() {
    kontrola_kroku(std::pow(10.,-2), "trapezy_1_TOL_2.txt", metTrapezow);
    kontrola_kroku(std::pow(10.,-5), "trapezy_1_TOL_5.txt", metTrapezow);

    kontrola_kroku(std::pow(10.,-2), "rk2_1_TOL_2.txt",metRK2);
    kontrola_kroku(std::pow(10.,-5),"rk2_1_TOL_5.txt",metRK2);

    return 0;
}

void kontrola_kroku(double TOL, std::string filename, std::function<std::array<double,2>(double,double,double,double)>metoda){
    double x0 = 0.01;
    double v0 = 0.;
    double dt0 = 1.;
    double dt = dt0;
    double S = 0.75;
    double p = 2.;
    double tmax = 40.;
    double alfa = 5.;

    double t = 0.;
    double xn = x0;
    double vn = v0;


    //double TOL = std::pow(10.,-2);
    FILE *fp1 = fopen(filename.c_str(),"w");

    do{
        //stawiamy dwa kroki dt
        double x_pierwszy = metoda(xn,vn,dt,alfa)[0];
        double v_pierwszy = metoda(xn,vn,dt,alfa)[1];

        double x_drugi = metoda(x_pierwszy,v_pierwszy,dt,alfa)[0];
        double v_drugi = metoda(x_pierwszy,v_pierwszy,dt,alfa)[1];

        //stawiamy jeden krok 2dt
        double x_2dt = metoda(xn,vn,2*dt,alfa)[0];
        double v_2dt = metoda(xn,vn,2*dt,alfa)[1];

        double Ex = (x_drugi - x_2dt)/(2*p - 1);
        double Ev = (v_drugi - v_2dt)/(2*p - 1);

        if(std::max(std::abs(Ex),std::abs(Ev)) < TOL){
            t = t + 2*dt;
            xn = x_drugi;
            vn = v_drugi;
            fprintf(fp1,"%12.6g\t%12.6g\t%12.6g\t%12.6g\n",t,dt,xn,vn);

//            std::cout<<t<<"\t"<<dt<<"\t"<<xn<<"\t"<<vn<<std::endl;
        }

        dt = dt * std::pow((S*TOL)/std::max(std::abs(Ex),std::abs(Ev)),1./(p+1.));

    }while(t < tmax);

}
