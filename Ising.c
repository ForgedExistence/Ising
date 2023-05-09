//Ising Model
/*Dónal McGrath ID:7012616
Prajit Baruah ID: 9097760
*/

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 2
#define N 100

double s[N][N]={0};
int n_particles = 0;
double acc = 0.5;
int mc_steps = 50000;
int measure = 1000; //How often we'll measure shit
int rand_i, rand_j;
double num;
double total=0;
int a, b;
double init_e = 0; 
double J = 1; //Ask if this needs to be changed. 
int T = 0.1; 
double mag = 0; 
double ct = 0.1; // temperature change delta 
double avg_energy = 0.0; 
double energy = 0; 
double avg_mag = 0; 
double e2 = 0; //energy square for specific heat 
double heat = 0; 
double beta = 0; 

/*
If initial code works, try and implement some of the improvements mentioned in the book. 
*/


int IsingPos(void){
    char buffer[128];
    sprintf(buffer, "State_step_Init.dat");
    //FILE* fp = fopen(buffer, "w");
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            s[i][j] = 2*(dsfmt_genrand()-0.5);
            if(s[i][j]<0){s[i][j] = -1.;}
            else{s[i][j] = 1.;}
            //fprintf(fp, "%d %d %lf\n", i, j, s[i][j]);
        }
    }
    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            for (int x = i-1; x<= i+1; x++){
                for(int y = j-1; y <= j+1; y++){
                    
                    if(x==i && y==j){continue;}
                    a=x;
                    b=y;
                    if(x==-1){a=N-1; }
                    if(x==N){a = 0;}
                    if(y==-1){b=N-1;} //Check if simpler way to enforce PBC
                    if(y==N){b=0;}

                    init_e += J*s[i][j]*s[a][b]; //Calculate initial energy of the system 
                }
            }
        }
    }
    //fclose(fp);
}

int interact(int i, int j){
    double de = 0; 
    energy = init_e; 
    beta = 1/T; 
    //printf("pos: %d %d old spin: %lf\n", i, j, s[i][j]);
    for(int n = i-1; n <= i+1;n++){
        for(int m = j-1; m <= j+1;m++){
            if(n==i && m==j){continue;}
            a=n;
            b=m;
            if(n==-1){a=N-1; }
            if(n==N){a = 0;}
            if(m==-1){b=N-1;}
            if(m==N){b=0;}
            
            de += 2*J*s[i][j]*s[a][b]; //Energy change 
        }
    } 

    if(de < 0.0 || dsfmt_genrand() < exp(-beta*de)){

        s[i][j] *= -1;
        energy += de; //to track total energy 
        return 1; 
    }

    else{return 0;}
}

void write_data(double step){

    FILE* fp = fopen("Isingdata.dat", "w");
    fprintf(fp,"%d\t%lf\t%lf\t%lf\n", T,avg_mag, avg_energy, heat);
    fclose(fp);
}

int main(int argc, char* argv[]){
    int div = (mc_steps-10000)/measure; 
    printf("%d\n", div);
    dsfmt_seed(time(NULL));
    struct stat st = {0};

    if(stat("data_Ising", &st)==-1){
       // mkdir("data_Ising");
        //printf("created folder\"data_Ising\"\n");
    }  
    IsingPos();

    while(T <= 5){

    for(int steps = 0; steps< mc_steps; steps++){
        rand_i = (int) (N*dsfmt_genrand());
        rand_j = (int) (N*dsfmt_genrand());

        //printf("pos: %d %d old spin: %lf\n", rand_i, rand_j, s[rand_i][rand_j]);
        interact(rand_i, rand_j);

        if(steps >= 10000 && steps % measure == 0){ //Confirm if system is actually in equillibrium

            avg_energy += energy;
            e2 += energy*energy;
            for(int i = 0; i<N; i++){
                for (int j = 0; j<N; j++){
                    avg_mag += s[i][j];
                }
            }
        }
    }
    avg_energy /= div; 
    avg_mag /= div;
    e2 /= div;
    heat = beta*beta*(e2 - (avg_energy*avg_energy))/(N*N);
    T += ct;
    }
}