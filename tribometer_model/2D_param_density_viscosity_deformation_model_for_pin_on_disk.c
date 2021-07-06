#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define PI acos(-1)
#define R 7.5e-3 // radius of the pin
#define omage 1.0e-3 // learning rate of iteration

#define omega 0.01 // iteration rate
#define M 100 //maximum node in x-axies
#define N 100 //maximum node in y-axies
#define R 7.5e-3 // pin's radius
#define PI acos(-1) //pi
#define Err_Pmin 1.0e-4 //convergence criteria
double uA = 0.0209,uB = 0.0; // velocity 100rpm [m/s] in x-axies
double u_avg = 0.0209/2;
double v = 0.0; // velocity in y-axies
double U = 1.0;
double b = 1.5e-3/2.0; // x-axies
double a = 1.5e-3/2.0; // y-axies
double k = 1.0; // b/a
double rho_0 = 1.82e3; // [kg/m3]
double eta_0 = 0.055; //[Pa.s]
double alpha_0 = 50.0e-9; // [Pa-1]
double hmin = 6e-9; //clearance [m]
double P_h = 1.01e5;
double DX = 2.0/M, DY = 2.0/N;
double beta_x = 0.5;
double Psi;
double E1=150.0e9,E2=110.0e9,E_dash; //Young's molecule [Pa]
FILE *fp1,*fp2,*fp3;

struct Tribo{
    double dP[N+1][M+1];
    double P[N+1][M+1];
    double H[N+1][M+1];
    double HS[N+1][M+1];
    double H0[N+1][M+1];
    double Delta[N+1][M+1];
    double Rho_bar[N+1][M+1];
    double Eta_bar[N+1][M+1];
    double X[M+1];
    double Y[N+1];
};

//calculate F(i,j)
double cal_F_ij(int i, int j, double Eta_bar[N+1][M+1], double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double Ax,Ay,Bx,By,C;
    double rho_H3_eta_1,rho_H3_eta_2,rho_H3_eta_3,rho_H3_eta_4,rho_H3_eta_5;
    double rho_H_U_1, rho_H_U_2, rho_H_U_3;
    double U = (uA+uB)/2/u_avg;
    double F_i_j;
    //Rho_H3_Eta{(i-1,j), (i+1,j), (i,j), (i,j-1), (i,j+1)}
    rho_H3_eta_1 = Rho_bar[j][i-1]*pow(H[j][i-1], 3)/Eta_bar[j][i-1];
    rho_H3_eta_2 = Rho_bar[j][i+1]*pow(H[j][i+1], 3)/Eta_bar[j][i+1];
    rho_H3_eta_3 = Rho_bar[j][i]*pow(H[j][i], 3)/Eta_bar[j][i];
    rho_H3_eta_4 = Rho_bar[j-1][i]*pow(H[j-1][i], 3)/Eta_bar[j-1][i];
    rho_H3_eta_5 = Rho_bar[j+1][i]*pow(H[j+1][i], 3)/Eta_bar[j+1][i];
    
    //RhoHU{(i-1,j), (i,j), (i+1,j)}
    rho_H_U_1 = Rho_bar[j][i-1]*H[j][i-1]*U;
    rho_H_U_2 = Rho_bar[j][i]*H[j][i]*U;
    rho_H_U_3 = Rho_bar[j][i+1]*H[j][i+1]*U;
    //Ax,Ay,Bx,By,C;
    Ax = P[j][i-1]*(rho_H3_eta_1+rho_H3_eta_3) - P[j][i]*(rho_H3_eta_1+rho_H3_eta_2+2*rho_H3_eta_3) + P[j][i+1]*(rho_H3_eta_2+rho_H3_eta_3);
    Ay = P[j-1][i]*(rho_H3_eta_4+rho_H3_eta_3) - P[j][i]*(rho_H3_eta_4+rho_H3_eta_5+2*rho_H3_eta_3) + P[j+1][i]*(rho_H3_eta_5+rho_H3_eta_3);
    Bx = (1-beta_x)*(rho_H_U_3-rho_H_U_2) + beta_x*(rho_H_U_2-rho_H_U_1);
    By = 0.0;
    C = 0.0;
    //F = 1/DX2 * Ax + k2 * 1/DY2 *Ay - Psi{Bx+kBy+C}
    F_i_j = 1/(2*DX*DX) * Ax + 1/(2*DY*DY) * Ay - Psi*(Bx+By+C);
    return F_i_j;
}
//Ro,Et,Dm
double cal_Ro_ij_kl(int i, int j, int k, int l, double P_ij){
    if (k==i && l==j) {
        return 0.6e-9*P_h/pow(1+1.7e-9*P_h*P_ij,2);
    }else{
        return 0.0;
    }
}
double cal_Et_ij_kl(int i, int j, int k, int l, double P_ij, double Eta_bar_ij){
    double z = alpha_0/(5.1e-9*(log(eta_0)+9.67));
    if (k == i && l == j) {
        return (log(eta_0)+9.67)*z*P_h*Eta_bar_ij/1.98e8*pow(1+P_ij*P_h/1.98e8, z-1);
    }else{
        return 0.0;
    }
}
double cal_D_mn(int m,int n){
    double dx,dy;
    double x_bar,y_bar,a_bar,b_bar;
    double xm,xp,ym,yp;
    double A,B,C,D;
    dx = b*DX;
    dy = a*DY;
    x_bar = m*dx;
    y_bar = n*dy;
    a_bar = dy/2.0;
    b_bar = dx/2.0;
    xm = x_bar-b_bar;
    xp = x_bar+b_bar;
    ym = y_bar-a_bar;
    yp = y_bar+a_bar;
    A = ym*log((xm+sqrt(pow(ym, 2)+pow(xm, 2)))/(xp+sqrt(pow(ym, 2)+pow(xp, 2))));
    B = yp*log((xp+sqrt(pow(yp, 2)+pow(xp, 2)))/(xm+sqrt(pow(yp, 2)+pow(xm, 2))));
    C = xp*log((yp+sqrt(pow(yp, 2)+pow(xp, 2)))/(ym+sqrt(pow(ym, 2)+pow(xp, 2))));
    D = xm*log((ym+sqrt(pow(ym, 2)+pow(xm, 2)))/(yp+sqrt(pow(yp, 2)+pow(xm, 2))));
    return A+B+C+D;
}
//delta(k,l)
double cal_delta_kl(int k, int l, double P[N+1][M+1]){
    double delta_kl = 0.0;
    int m,n;
    double D_mn;
    for (int j=0; j<=N; j++) {
        for (int i=0; i<=M; i++) {
            m = abs(k-i);
            n = abs(l-j);
            D_mn = cal_D_mn(m, n);
            delta_kl+= (2*P_h/(PI*E_dash)) * P[j][i] * D_mn;
        }
    }
    return delta_kl;
}
// M(ij,kl) , N(ij,kl)
double cal_M_ij_kl(int i, int j, int k, int l, double Eta_bar[N+1][M+1], double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double H3_Eta_ij, Rho_H3_Eta2_ij, Rho_H2_Eta_ij;
    double Et_ij_kl, Rh_ij_kl, D_mn;
    int m,n;
    m = abs(k-i);
    n = abs(l-j);
    Rh_ij_kl = cal_Ro_ij_kl(i, j, k, l, P[j][i]);
    Et_ij_kl = cal_Et_ij_kl(i, j, k, l, P[j][i], Eta_bar[j][i]);
    D_mn = cal_D_mn(m, n);
    H3_Eta_ij = pow(H[j][i], 3)/Eta_bar[j][i];
    Rho_H2_Eta_ij = Rho_bar[j][i]*pow(H[j][i], 2)/Eta_bar[j][i];
    Rho_H3_Eta2_ij = Rho_bar[j][i]*pow(H[j][i], 3)/pow(Eta_bar[j][i], 2);
    
    return H3_Eta_ij*Rh_ij_kl - Et_ij_kl*Rho_H3_Eta2_ij + 3*Rho_H2_Eta_ij * D_mn;
}
double cal_Nx_ij_kl(int i, int j, int k, int l, double Rho_bar[N+1][M+1],double H[N+1][M+1], double P[N+1][M+1]){
    double Ro_ij_kl, D_mn;
    double HU_ij, Rho_U_ij;
    int m,n;
    m = abs(k-i);
    n = abs(l-j);
    HU_ij = H[j][i]*U;
    Rho_U_ij = Rho_bar[j][i]*U;
    D_mn = cal_D_mn(m, n);
    Ro_ij_kl = cal_Ro_ij_kl(i,j, k, l, P[j][i]);
    
    return HU_ij*Ro_ij_kl + Rho_U_ij*D_mn;
}

//J(ij,ij) = dF(i,j)/dP(i,j)
double cal_J_ij_ij(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double rho_H3_eta_im1j,rho_H3_eta_ij,rho_H3_eta_ip1j,rho_H3_eta_ijp1,rho_H3_eta_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j, i, j, Eta_bar, Rho_bar, H, P);
    M2 = cal_M_ij_kl(i+1, j, i, j, Eta_bar, Rho_bar, H, P);
    M3 = cal_M_ij_kl(i, j, i, j, Eta_bar, Rho_bar, H, P);
    M4 = cal_M_ij_kl(i, j-1, i, j, Eta_bar, Rho_bar, H, P);
    M5 = cal_M_ij_kl(i, j+1, i, j, Eta_bar, Rho_bar, H, P);
    Nx1 = cal_Nx_ij_kl(i-1, j, i, j, Rho_bar, H, P);
    Nx2 = cal_Nx_ij_kl(i, j, i, j, Rho_bar, H, P);
    Nx3 = cal_Nx_ij_kl(i+1, j, i, j, Rho_bar, H, P);
    
    rho_H3_eta_ij = Rho_bar[j][i]*pow(H[j][i], 3)/Eta_bar[j][i];
    rho_H3_eta_im1j = Rho_bar[j][i-1]*pow(H[j][i-1], 3)/Eta_bar[j][i-1];
    rho_H3_eta_ip1j = Rho_bar[j][i+1]*pow(H[j][i+1], 3)/Eta_bar[j][i+1];
    rho_H3_eta_ijm1 = Rho_bar[j-1][i]*pow(H[j-1][i], 3)/Eta_bar[j-1][i];
    rho_H3_eta_ijp1 = Rho_bar[j+1][i]*pow(H[j+1][i], 3)/Eta_bar[j+1][i];
    
    A1 = 1/(2*DX*DX) * (P[j][i-1]*(M3+M1) - P[j][i]*(M2+M3*2+M1) + P[j][i+1]*(M2+M3) - (rho_H3_eta_im1j+rho_H3_eta_ip1j+2*rho_H3_eta_ij));
    A2 = k*k/(2*DY*DY) * (P[j-1][i]*(M4+M3) - P[j][i]*(M4+M3*2+M5) + P[j+1][i] * (M5+M3) - (rho_H3_eta_ijm1+rho_H3_eta_ijp1+2*rho_H3_eta_ij));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C =0.0;
    
    return A1+A2-B1-B2-C;
}
//J(ij,i-1j) = dF(i,j)/dP(i-1,j)
double cal_J_ij_im1j(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double rho_H3_eta_im1j,rho_H3_eta_ij,rho_H3_eta_ip1j,rho_H3_eta_ijp1,rho_H3_eta_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j, i-1, j, Eta_bar, Rho_bar, H, P);
    M2 = cal_M_ij_kl(i+1, j, i-1, j, Eta_bar, Rho_bar, H, P);
    M3 = cal_M_ij_kl(i, j, i-1, j, Eta_bar, Rho_bar, H, P);
    M4 = cal_M_ij_kl(i, j-1, i-1, j, Eta_bar, Rho_bar, H, P);
    M5 = cal_M_ij_kl(i, j+1, i-1, j, Eta_bar, Rho_bar, H, P);
    Nx1 = cal_Nx_ij_kl(i-1, j, i-1, j, Rho_bar, H, P);
    Nx2 = cal_Nx_ij_kl(i, j, i-1, j, Rho_bar, H, P);
    Nx3 = cal_Nx_ij_kl(i+1, j, i-1, j, Rho_bar, H, P);
    
    rho_H3_eta_ij = Rho_bar[j][i]*pow(H[j][i], 3)/Eta_bar[j][i];
    rho_H3_eta_im1j = Rho_bar[j][i-1]*pow(H[j][i-1], 3)/Eta_bar[j][i-1];
    rho_H3_eta_ip1j = Rho_bar[j][i+1]*pow(H[j][i+1], 3)/Eta_bar[j][i+1];
    rho_H3_eta_ijm1 = Rho_bar[j-1][i]*pow(H[j-1][i], 3)/Eta_bar[j-1][i];
    rho_H3_eta_ijp1 = Rho_bar[j+1][i]*pow(H[j+1][i], 3)/Eta_bar[j+1][i];
    
    A1 = 1/(2*DX*DX) * (P[j][i-1]*(M3+M1) - P[j][i]*(M2+M3*2+M1) + P[j][i+1]*(M2+M3) + (rho_H3_eta_im1j+rho_H3_eta_ij));
    A2 = k*k/(2*DY*DY) * (P[j-1][i]*(M4+M3) - P[j][i]*(M4+M3*2+M5) + P[j+1][i] * (M5+M3));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C =0.0;
    
    return A1+A2-B1-B2-C;
}
//J(ij,i+1j) = dF(i,j)/dP(i+1,j)
double cal_J_ij_ip1j(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double rho_H3_eta_im1j,rho_H3_eta_ij,rho_H3_eta_ip1j,rho_H3_eta_ijp1,rho_H3_eta_ijm1;
    double A1,A2,B1,B2,C;
    M1 = cal_M_ij_kl(i-1, j, i+1, j, Eta_bar, Rho_bar, H, P);
    M2 = cal_M_ij_kl(i+1, j, i+1, j, Eta_bar, Rho_bar, H, P);
    M3 = cal_M_ij_kl(i, j, i+1, j, Eta_bar, Rho_bar, H, P);
    M4 = cal_M_ij_kl(i, j-1, i+1, j, Eta_bar, Rho_bar, H, P);
    M5 = cal_M_ij_kl(i, j+1, i+1, j, Eta_bar, Rho_bar, H, P);
    Nx1 = cal_Nx_ij_kl(i-1, j, i+1, j, Rho_bar, H, P);
    Nx2 = cal_Nx_ij_kl(i, j, i+1, j, Rho_bar, H, P);
    Nx3 = cal_Nx_ij_kl(i+1, j, i+1, j, Rho_bar, H, P);
    
    rho_H3_eta_ij = Rho_bar[j][i]*pow(H[j][i], 3)/Eta_bar[j][i];
    rho_H3_eta_im1j = Rho_bar[j][i-1]*pow(H[j][i-1], 3)/Eta_bar[j][i-1];
    rho_H3_eta_ip1j = Rho_bar[j][i+1]*pow(H[j][i+1], 3)/Eta_bar[j][i+1];
    rho_H3_eta_ijm1 = Rho_bar[j-1][i]*pow(H[j-1][i], 3)/Eta_bar[j-1][i];
    rho_H3_eta_ijp1 = Rho_bar[j+1][i]*pow(H[j+1][i], 3)/Eta_bar[j+1][i];
    
    A1 = 1/(2*DX*DX) * (P[j][i-1]*(M3+M1) - P[j][i]*(M2+M3*2+M1) + P[j][i+1]*(M2+M3) + (rho_H3_eta_ip1j+rho_H3_eta_ij));
    A2 = k*k/(2*DY*DY) * (P[j-1][i]*(M4+M3) - P[j][i]*(M4+M3*2+M5) + P[j+1][i] * (M5+M3));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C =0.0;
    return A1+A2-B1-B2-C;
}
//J(ij,ij-1) = dF(i,j)/dP(i,j-1)
double cal_J_ij_ijm1(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double rho_H3_eta_im1j,rho_H3_eta_ij,rho_H3_eta_ip1j,rho_H3_eta_ijp1,rho_H3_eta_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j, i, j-1, Eta_bar, Rho_bar, H, P);
    M2 = cal_M_ij_kl(i+1, j, i, j-1, Eta_bar, Rho_bar, H, P);
    M3 = cal_M_ij_kl(i, j, i, j-1, Eta_bar, Rho_bar, H, P);
    M4 = cal_M_ij_kl(i, j-1, i, j-1, Eta_bar, Rho_bar, H, P);
    M5 = cal_M_ij_kl(i, j+1, i, j-1, Eta_bar, Rho_bar, H, P);
    Nx1 = cal_Nx_ij_kl(i-1, j, i, j-1, Rho_bar, H, P);
    Nx2 = cal_Nx_ij_kl(i, j, i, j-1, Rho_bar, H, P);
    Nx3 = cal_Nx_ij_kl(i+1, j, i, j-1, Rho_bar, H, P);
    
    rho_H3_eta_ij = Rho_bar[j][i]*pow(H[j][i], 3)/Eta_bar[j][i];
    rho_H3_eta_im1j = Rho_bar[j][i-1]*pow(H[j][i-1], 3)/Eta_bar[j][i-1];
    rho_H3_eta_ip1j = Rho_bar[j][i+1]*pow(H[j][i+1], 3)/Eta_bar[j][i+1];
    rho_H3_eta_ijm1 = Rho_bar[j-1][i]*pow(H[j-1][i], 3)/Eta_bar[j-1][i];
    rho_H3_eta_ijp1 = Rho_bar[j+1][i]*pow(H[j+1][i], 3)/Eta_bar[j+1][i];
    
    A1 = 1/(2*DX*DX) * (P[j][i-1]*(M3+M1) - P[j][i]*(M2+M3*2+M1) + P[j][i+1]*(M2+M3));
    A2 = k*k/(2*DY*DY) * (P[j-1][i]*(M4+M3) - P[j][i]*(M4+M3*2+M5) + P[j+1][i]*(M5+M3)+(rho_H3_eta_ijm1+rho_H3_eta_ij));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C =0.0;
    
    return A1+A2-B1-B2-C;
}
//J(ij,ij+1) = dF(i,j)/dP(i,j+1)
double cal_J_ij_ijp1(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double rho_H3_eta_im1j,rho_H3_eta_ij,rho_H3_eta_ip1j,rho_H3_eta_ijp1,rho_H3_eta_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j, i, j+1, Eta_bar, Rho_bar, H, P);
    M2 = cal_M_ij_kl(i+1, j, i, j+1, Eta_bar, Rho_bar, H, P);
    M3 = cal_M_ij_kl(i, j, i, j+1, Eta_bar, Rho_bar, H, P);
    M4 = cal_M_ij_kl(i, j-1, i, j+1, Eta_bar, Rho_bar, H, P);
    M5 = cal_M_ij_kl(i, j+1, i, j+1, Eta_bar, Rho_bar, H, P);
    Nx1 = cal_Nx_ij_kl(i-1, j, i, j+1, Rho_bar, H, P);
    Nx2 = cal_Nx_ij_kl(i, j, i, j+1, Rho_bar, H, P);
    Nx3 = cal_Nx_ij_kl(i+1, j, i, j+1, Rho_bar, H, P);
    
    rho_H3_eta_ij = Rho_bar[j][i]*pow(H[j][i], 3)/Eta_bar[j][i];
    rho_H3_eta_im1j = Rho_bar[j][i-1]*pow(H[j][i-1], 3)/Eta_bar[j][i-1];
    rho_H3_eta_ip1j = Rho_bar[j][i+1]*pow(H[j][i+1], 3)/Eta_bar[j][i+1];
    rho_H3_eta_ijm1 = Rho_bar[j-1][i]*pow(H[j-1][i], 3)/Eta_bar[j-1][i];
    rho_H3_eta_ijp1 = Rho_bar[j+1][i]*pow(H[j+1][i], 3)/Eta_bar[j+1][i];
    
    A1 = 1/(2*DX*DX) * (P[j][i-1]*(M3+M1) - P[j][i]*(M2+M3*2+M1) + P[j][i+1]*(M2+M3));
    A2 = k*k/(2*DY*DY) * (P[j-1][i]*(M4+M3) - P[j][i]*(M4+M3*2+M5) + P[j+1][i]*(M5+M3) + (rho_H3_eta_ijp1+rho_H3_eta_ij));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C =0.0;
    
    return A1+A2-B1-B2-C;
}

void init_param(struct Tribo *node){
    double h0,hs,delta;
    double p,P;
    double x,y,x0=-b,y0=-a,dx,dy;
    double z;
    for (int i=0; i<=M; i++) {
        node->X[i] = -1.0+i*DX;
    }
    for (int j=0; j<=N; j++) {
        node->Y[j] = -1.0+j*DY;
    }
    for (int j=0; j<=N; j++) {
        for (int i=0; i<=M; i++) {
            //pressure
            p = 0.0;
            P = p/P_h;
            node->P[j][i] = p/P_h;
            node->dP[j][i] = 0.0;
            //film thickness
            dx = DX*b;
            dy = DY*a;
            x = x0+i*dx;
            y = y0+j*dy;
            h0 = hmin;
            hs = R-sqrt(pow(R, 2)-pow(x, 2)-pow(y, 2));
            delta = 0.0;
            node->HS[j][i] = hs/(b*b/R);
            node->H0[j][i] = h0/(b*b/R);
            node->Delta[j][i] = delta/(b*b/R);
            node->H[j][i] = (h0+hs+delta)/(b*b/R);
            //density
            node->Rho_bar[j][i] = 1+(0.6e-9*P_h*P)/(1+1.7e-9*P_h*P);
            //viscosity
            z = alpha_0/(5.1e-9*(log(eta_0)+9.67));
            node->Eta_bar[j][i] = exp((log(eta_0)+9.67) * (-1 + pow(1+P_h*P/1.98e8, z)));
        }
    }
    
}
void print_matrix(double Array[N+1][M+1]){
    for (int j=0; j<=N; j++) {
        for (int i=0; i<=M; i++) {
            printf(" %.5lf",Array[j][i]);
        }printf("\n");
    }printf("\n");
}
void update_param(struct Tribo *node){
    double F,J1,J2,J3,J4,J5;
    double delta_kl;
    double z;
    E_dash = 1/(1/E1+1/E2);
    z = alpha_0/(5.1e-9*(log(eta_0)+9.67));
    for (int j=1; j < N ; j++) {
        for (int i = 1; i<M; i++) {
            F = cal_F_ij(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P); //F(i,j)
            J1 = cal_J_ij_im1j(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P); //J(i-1,j)
            J2 = cal_J_ij_ip1j(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P); //J(i+1,j)
            J3 = cal_J_ij_ij(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P); //J(i,j)
            J4 = cal_J_ij_ijm1(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P); //J(i,j-1)
            J5 = cal_J_ij_ijp1(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P); //J(i,j+1)
            node->dP[j][i] = -(F + J1*node->dP[j][i-1] + J2*node->dP[j][i+1] + J4*node->dP[j-1][i] + J5*node->dP[j+1][i]) / J3;
            node->P[j][i] = node->P[j][i] + omega*node->dP[j][i];
            if (node->P[j][i] < 0.0) {
                node->P[j][i] = 0.0;
            }
            //update viscosity
            node->Eta_bar[j][i] = exp((log(eta_0)+9.67) * (-1 + pow(1+P_h*node->P[j][i]/1.98e8, z)));
            //update density
            node->Rho_bar[j][i] = 1+(0.6e-9*P_h*node->P[j][i])/(1+1.7e-9*P_h*node->P[j][i]);
            //update film thickness
            for (int l= 0; l <= N; l++) {
                for (int k = 0; k <= M; k++ ) {
                    delta_kl = cal_delta_kl(k, l, node->P);
                    node->Delta[l][k] = delta_kl/(b*b/R);
                    node->H[l][k] = node->H0[l][k]+node->HS[l][k]+node->Delta[l][k];
                }
            }
        }
    }
}
//copy array
void copy_array(double new_array[N+1][M+1],double array[N+1][M+1]){
    for (int j=0; j<=N; j++) {
        for (int i=0; i<=M; i++) {
            new_array[j][i] = array [j][i];
        }
    }
}
int main(){
    struct Tribo node;
    double P_pervious[N+1][M+1], P_new[N+1][M+1];
    double Error_P;
    int loop=0;
    Psi = 12*u_avg*eta_0*R*R/pow(b, 3)/P_h;
    init_param(&node);
    print_matrix(node.HS);
    double diff_P = 0.0, sum_P = 0.0;
    double pMax = 0.0, hMin = 1e10;
    //loop
    do {
        loop++;
        printf("the loop is : %d-th\n",loop);
        copy_array(P_pervious, node.P);
        update_param(&node);

        copy_array(P_new, node.P);

        for (int j=0; j<=N; j++) {
            for (int i=0; i<=M; i++) {
                if (node.P[j][i]>0.0) {
                    sum_P+=fabs(P_new[j][i]);
                    diff_P+=fabs(P_new[j][i]-P_pervious[j][i]);
                }
                if(node.P[j][i]>pMax) pMax = node.P[j][i];
                if(node.H[j][i]<hMin) hMin = node.H[j][i];
                printf("%15.10lf ",node.Delta[j][i]);
            }printf("\n");
        }printf("\n");
        Error_P = diff_P/sum_P;
        printf("the error of pressure is : %.5lf\n",Error_P);

    } while (Error_P>Err_Pmin || loop < 100);
    
    fp1 = fopen("/Users/liyurun/Desktop/simulationdata/2D_Parameter_Rho_Eta_deform_Pressure.csv", "w+");
    for (int j = 0; j<N+1; j++){
        for (int i = 0; i< M+1; i++){
            fprintf(fp1, "%.5lf%s",node.P[j][i]*P_h,
                    (i<M+1-1?",":""));
        }
        fprintf(fp1,"\n");
    }
    fclose(fp1);

    fp2 = fopen("/Users/liyurun/Desktop/simulationdata/2D_Parameter_Rho_Eta__deform_FilmThickness.csv", "w+");
    for (int j = 0; j<N+1; j++){
        for (int i = 0; i< M+1; i++){
            fprintf(fp2, "%.10lf%s",node.H[j][i]*(b*b/R),
                    (i<M+1-1?",":""));
        }
        fprintf(fp2,"\n");
    }
    fclose(fp2);
    
    fp3 = fopen("/Users/liyurun/Desktop/simulationdata/2D_Parameter_Rho_Eta_deform_Deformation.csv", "w+");
    for (int j = 0; j<N+1; j++){
        for (int i = 0; i< M+1; i++){
            fprintf(fp3, "%.10lf%s",node.Delta[j][i]*(b*b/R),
                    (i<M+1-1?",":""));
        }
        fprintf(fp3,"\n");
    }
    fclose(fp3);
}



