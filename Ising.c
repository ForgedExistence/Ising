//Ising Model
/*Dónal McGrath ID:7012616
Prajit Baruah*/

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
#define N 8

double s[N][N]={0};
int n_particles = 0;
double acc = 0.5;
int mc_steps = 100000;
int rand_i, rand_j;
double num;
double total=0;
int a, b;


int IsingPos(void){
    char buffer[128];
    sprintf(buffer, "State_step_Init.dat");
    FILE* fp = fopen(buffer, "w");
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            s[i][j] = 2*(dsfmt_genrand()-0.5);
            if(s[i][j]<0){s[i][j] = -1.;}
            else{s[i][j] = 1.;}
            fprintf(fp, "%d %d %lf\n", i, j, s[i][j]);
        }
    }
    fclose(fp);
}

int interact(int i, int j){
    num = 0.;
    printf("pos: %d %d old spin: %lf\n", i, j, s[i][j]);
    for(int n = i-1; n <= i+1;n++){
        for(int m = j-1; m <= j+1;m++){
            if(n==i && m==j){continue;}
            a=n;
            b=m;
            if(n==-1){a=N-1; }
            if(n==N){a = 0;}
            if(m==-1){b=N-1;}
            if(m==N){b=0;}
            num += s[a][b];
            //printf("%d %d\n%d %d\n%lf\n", a, b, n, m, s[a][b]);
        }
    } 
    if(num <0){s[i][j] = -1.;}
    if(num >0){s[i][j] = 1.;}
    //printf("num: %lf\n", num);
    //printf("pos: %d %d updated spin: %lf\n", i, j, s[i][j]);
}

void write_data(double step){
    char buffer[128];
    sprintf(buffer, "Ising/State_step%07lf.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < N; ++d) fprintf(fp, "%d %d", n, d);
        fprintf(fp, "%lf\n", s[n][d]);
    }
    fclose(fp);
}

int main(int argc, char* argv[]){
    dsfmt_seed(time(NULL));
    struct stat st = {0};

    if(stat("data_Ising", &st)==-1){
        mkdir("data_Ising");
        printf("created folder\"data_Ising\"\n");
    }  
    IsingPos();
    for(int steps = 0; steps< mc_steps; steps++){
        rand_i = (int) (N*dsfmt_genrand());
        rand_j = (int) (N*dsfmt_genrand());
        //printf("pos: %d %d old spin: %lf\n", rand_i, rand_j, s[rand_i][rand_j]);
        interact(rand_i, rand_j);
    }
    for(int i =0; i<N; i++){
        for(int j=0; j<N; j++){
            //printf("%d %d %lf\n", i, j, s[i][j]);
            total += s[i][j]; 
        }
    }
    printf("total mag: %lf", total/(N*N));



}