#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define PI acos(-1)
#define Err_Pmin 1.0e-3
#define R 7.5e-3
#define omage 1.0e-3

double uA = 0.0209; //velocity
double u_avg =0.0209/2.0;
double L = 1.5e-3,b = 1.5e-3/2.0;
double rho = 1.82e3,eta = 0.055;
double hmin = 6.0e-9;
double P_H = 1.01e5;
double DX = 2.0/N;
double Psi,beta;
double beta_bar;

struct Tribo{
    double dP[N+1];
    double P[N+1];
    double H[N+1];
    double HS[N+1];
    double H0[N+1];
    double Delta[N+1];
    double Rho[N+1];
    double Eta[N+1];
};

void initialization(struct Tribo *node){
    double h0=hmin,hs,delta,h,dx,x;
    double p;
    dx = L/(N*1.0);
    for (int i=0; i<=N; i++) {
        node->P[i] = 0.0/P_H;
        node->dP[i] = 0.0;
        node->Eta[i] = 1.0;
        node->H0[i] = hmin/(b*b/R);
        x = i*dx;
        hs = R-sqrt(pow(R, 2)-pow(x-b, 2));
        node->HS[i] = hs/(b*b/R);
        node->Delta[i]=0.0;
        node->H[i] = node->HS[i]+node->H0[i]+node->Delta[i];
        node->Rho[i] = 1.0;
        node->Eta[i] = 1.0;
        
    }
}
double cal_F_i(int i, double H[N+1],double P[N+1]){
    double A,C;
    double H3_1,H3_2,H3_3;
    H3_1 = pow(H[i-1], 3);
    H3_2 = pow(H[i], 3);
    H3_3 = pow(H[i+1], 3);
    
    A = (H3_3+H3_2)*P[i+1] - (H3_3+2*H3_2+H3_1)*P[i] + (H3_2-H3_1)*P[i-1];
    C = (1-beta)*(H[i+1]-H[i]) + beta*(H[i]-H[i-1]);
    return 1/(2*DX*DX)*A - (Psi/DX)*C;
}

double cal_J_im1(int i, double H[N+1]){
    double H3_1,H3_2,H3_3;
    H3_1 = pow(H[i-1], 3);
    H3_2 = pow(H[i], 3);
    H3_3 = pow(H[i+1], 3);
    
    return 1/(2*DX*DX)*(H3_1+H3_2);
    
}
double cal_J_i(int i, double H[N+1]){
    double H3_1,H3_2,H3_3;
    H3_1 = pow(H[i-1], 3);
    H3_2 = pow(H[i], 3);
    H3_3 = pow(H[i+1], 3);
    
    return -1/(2*DX*DX)*(H3_1+H3_2*2+H3_3);
    
}
double cal_J_ip1(int i, double H[N+1]){
    double H3_1,H3_2,H3_3;
    H3_1 = pow(H[i-1], 3);
    H3_2 = pow(H[i], 3);
    H3_3 = pow(H[i+1], 3);
    
    return 1/(2*DX*DX)*(H3_2+H3_3);
    
}

void CopyArray(double Array_old[N+1],double Array_new[N+1]){
    for (int i=0; i<=N; i++) {
        Array_new[i] = Array_old[i];
    }
}

void updateNode(struct Tribo *node){
    double F,J1,J2,J3,dP;
    for (int i=1; i<N; i++) {
        F = cal_F_i(i, node->H, node->P);
        J1 = cal_J_im1(i, node->H);
        J2 = cal_J_i(i, node->H);
        J3 = cal_J_ip1(i, node->H);
        
        dP = (-F-J1*node->dP[i-1]-J3*node->dP[i+1])/J2;
        node->P[i] = node->P[i] + omage*dP;
        if (node->P[i]<0.0) {
            node->P[i] = 0.0;
        }
    }
}

int main() {
    struct Tribo node;
    double P_pervious[N+1],P_new[N+1];
    double Error_P;
    beta = 0.5;
    Psi = 12.0*u_avg*eta*pow(R, 2)/(P_H*pow(b, 3));
    initialization(&node);
    double diff_P=0.0,sum_P=0.0;
    double pMax=0.0,hMin=10000.0;
    //loop
    int loop=0;
    do{
        loop++;
        printf("the loop is : No.%d\nThe pressure is : ",loop);
        CopyArray(node.P, P_pervious);
        updateNode(&node);
        CopyArray(node.P, P_new);
        
        for (int i = 0; i<=N; i++) {
            if (node.P[i]>0.0) {
                sum_P+=fabs(P_new[i]);
                diff_P+=fabs(P_new[i]-P_pervious[i]);
            }
            if(P_new[i]>pMax) pMax = P_new[i];
            if(node.H[i]*(b*b/R) < hMin) hMin = node.H[i]*(b*b/R);
            printf("%.5lf ",node.P[i]*P_H);
        }printf("\n");
        Error_P = diff_P/sum_P;
        printf("The Error is %.10lf\n",Error_P);
    }while(Error_P > Err_Pmin || loop<=100);
    return 0;
}
