#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>


//// Zadanie 1
void schematEulera(double dt, std::string filename, bool metAnalityczna);

//// Zadanie 2
void metoda_rk2(double dt, std::string filename);

//// Zadanie 3
void metoda_rk4(double dt, std::string filename);

//// Zadanie 4
void RLC(double k, std::string filename);

//// fun pomocnicze
double funY(double t, double lambda){return std::exp(t * lambda);}
double funV(double omega, double t){return 10 * sin(omega * t);}



int main() {
    //Zadanie 1
    schematEulera(0.01,"euler_1", false);
    schematEulera(0.1,"euler_2", false);
    schematEulera(1.,"euler_3", false);
    schematEulera(0.01,"none", true);

    //Zadanie 2
    metoda_rk2(0.01, "rk2_1");
    metoda_rk2(0.1, "rk2_2");
    metoda_rk2(1., "rk2_3");

    //Zadanie 3
    metoda_rk4(0.01,"rk4_1");
    metoda_rk4(0.1, "rk4_2");
    metoda_rk4(1., "rk4_3");

    //Zadanie 4
    RLC(0.5, "rlc_0_5.txt");
    RLC(0.8, "rlc_0_8.txt");
    RLC(1., "rlc_1.txt");
    RLC(1.2,"rlc_1_2.txt");

    return 0;
}
//metAnalityczna = true - zwraca dane do wykresu dokladnego
void schematEulera(double dt, std::string filename, bool metAnalityczna){
    double lambda = -1.;
    double t0 = 0.;
    double t1 = 5.;
    int n = (t1 - t0)/dt;
    auto *y = new double [n+1];
    y[0]=1.;

    if(metAnalityczna){
        FILE *metAnalityczna = fopen("metAnalityczna.txt","w");
        for(int i=0; i<n; i++){
            fprintf(metAnalityczna,"%12.6g\t%12.6g\n",1.*i*dt,funY(1.*i*dt,lambda));
        }
        fclose(metAnalityczna);
    }else{

        std::string filename2 = filename + "_blad.txt";
        filename += ".txt";
        FILE *fp = fopen(filename.c_str(),"w");
        FILE *fp2 = fopen(filename2.c_str(),"w");
        for(int i=0; i<n; i++){
            double x = 1. * dt * i;
            double y_wynik = y[i] + dt * lambda * y[i];
            y[i+1] = y_wynik;
            //std::cout<<y[i]<<"\n";
            fprintf(fp,"%12.6g\t%12.6g\n",x,y[i]);
            fprintf(fp2,"%12.6g\t%12.6g\n",x,(y[i] - funY(x,lambda)));
        }
        if(filename == "euler_3.txt"){
            fprintf(fp,"%12.6g\t%12.6g\n",5*dt,y[5]);
            fprintf(fp2,"%12.6g\t%12.6g\n",5*dt,(y[5] - funY(5*dt,lambda)));

        }
        fclose(fp);
        fclose(fp2);
    }
    delete [] y;
}

void metoda_rk2(double dt, std::string filename){
    double lambda = -1.;
    double t0 = 0.;
    double t1 = 5.;
    int n = (t1 - t0)/dt;
    auto *y = new double [n+1];
    y[0]=1.;
    std::string filename2 = filename + "_blad.txt";
    filename += ".txt";

    FILE *fp = fopen(filename.c_str(),"w");
    FILE *fp2 = fopen(filename2.c_str(),"w");
    for(int i=0; i<n; i++){
        double x = 1. * i * dt;
        double k1 = lambda * y[i];
        double k2 = lambda * (y[i] + dt * k1);
        double y_wynik = y[i] + 1. * dt/2 * (k1 + k2);
        y[i+1]=y_wynik;

        fprintf(fp,"%12.6g\t%12.6g\n",x,y[i]);
        fprintf(fp2,"%12.6g\t%12.6g\n",x,y[i] - funY(x,lambda));
    }
    if(filename == "rk2_3.txt"){
        fprintf(fp,"%12.6g\t%12.6g\n",5*dt,y[5]);
        fprintf(fp2,"%12.6g\t%12.6g\n",5*dt,y[5] - funY(5*dt,lambda));
    }
    
    fclose(fp);
    fclose(fp2);
    delete [] y;
}

void metoda_rk4(double dt, std::string filename){
    double lambda = -1.;
    double t0 = 0.;
    double t1 = 5.;
    int n = (t1 - t0)/dt;
    auto *y = new double [n+1];
    y[0]=1.;
    std::string filename2 = filename + "_blad.txt";
    filename += ".txt";

    FILE *fp = fopen(filename.c_str(),"w");
    FILE *fp2 = fopen(filename2.c_str(),"w");
    for(int i=0 ;i<n; i++){
        double x = 1. * i * dt;
        double k1 = lambda * y[i];
        double k2 = lambda * (y[i] + dt/2 * k1);
        double k3 = lambda * ( y[i] + dt/2 * k2 );
        double k4 = lambda * ( y[i] + dt * k3);
        double y_wynik = y[i] + 1. * dt/6 * (k1 + 2 * k2 + 2 * k3 + k4);
        y[i+1]=y_wynik;

        fprintf(fp,"%12.6g\t%12.6g\n",x,y[i]);
        fprintf(fp2,"%12.6g\t%12.6g\n",x,y[i] - funY(x,lambda));
    }

    if(filename == "rk4_3.txt"){
        fprintf(fp,"%12.6g\t%12.6g\n",5*dt,y[5]);
        fprintf(fp2,"%12.6g\t%12.6g\n",5*dt,(y[5] - funY(5*dt,lambda)));
    }

    fclose(fp);
    fclose(fp2);
    delete [] y;
}

void RLC(double k, std::string filename){

    FILE *fp = fopen(filename.c_str(),"w");

    double dt = std::pow(10,-4);
    double R = 100;
    double L = 0.1;
    double C = 0.001;
    double omega = 1/std::pow(L*C,0.5);
    double To = 2 * M_PI / omega;
    int n = (4*To - 0.)/dt;

    auto * Q = new double [n+1];
    auto * I = new double [n+1];
    Q[0] = 0.;
    I[0] = 0.;

    double omegaT = k * omega;

    for(int i=0; i<n; i++){
        double x = dt * i;
        double k1_Q = I[i];
        double k1_I = funV(omegaT,x)/L - 1 / (L*C) * Q[i] - R/L * I[i];

        double k2_Q = I[i] + 1.*dt/2 * k1_I;
        double k2_I = funV(omegaT,(i + 0.5)*dt) / L - 1/(L*C) * (Q[i] + dt*0.5 * k1_Q) - R/L * (I[i] + dt*0.5 * k1_I);

        double k3_Q = I[i] + dt*0.5*k2_I;
        double k3_I = funV(omegaT,(i + 0.5)*dt)/L - 1/(L*C) * (Q[i] + dt*k3_Q) - R/L *(I[i] + dt*k3_I);

        double k4_Q = I[i] + dt * k3_I;
        double k4_I = funV(omegaT,(i+1)*dt)/L - 1/(L*C) * (Q[i] + dt*k3_Q) - R/L * (I[i] + dt*k3_I);

        Q[i+1] = Q[i] + dt/6 * (k1_Q + 2*k2_Q + 2*k3_Q + k4_Q);
        I[i+1] = I[i] + dt/6 * (k1_I + 2*k2_I + 2*k3_I + k4_I);

        fprintf(fp,"%12.6g\t%12.6g\t%12.6g\n",x, Q[i],I[i]);
    }

    delete [] I;
    delete [] Q;
    fclose(fp);
}
