//
//  main.c
//  Tribometer_rhe_theta
//
//  Created by YURUN LI on 2021/06/23.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI acos(-1)
#define omega 1.0e-3 // learning rate of iteration
#define M 10 //maximum node in x-axies
#define N 10 //maximum node in y-axies
#define PI acos(-1) //pi
#define Err_thetaMin 1.0e-4 //convergence criteria
double R = 7.5e-3; // pin's radius
double uA = 0.0209,uB = 0.0; // velocity 100rpm [m/s] in x-axies
double u_avg = 0.0209/2;
double v = 0.0; // velocity in y-axies
double U = 1.0;
double b = 20e-6; // x-axies
double a = 20e-6; // y-axies
double k = 1.0; // b/a
double rho_0 = 1.82e3; // [kg/m3]
double eta_0 = 0.055; //[Pa.s]
double alpha_0 = 50.0e-9; // [Pa-1]
double h_init = 1.7e-9; //initial film thickness
double h0 = 4e-9; // clearance [m]
double P_h = 1.01e5;
double p_c = 0.02e6; //cavitation pressure
double dx,dy,DX,DY; // unit distance
double beta_x = 0.5;
double Psi,beta_bar;
double E1=150.0e9,E2=110.0e9,E_dash; //Young's molecule [Pa]
double beta = 1.7e9; // bulk modulus of PFPE
FILE *fp1,*fp2,*fp3;

struct Tribo{
    double dP[N+1][M+1];
    double P[N+1][M+1];
    double H[N+1][M+1];
    double HS[N+1][M+1];
    double H0;
    double Delta[N+1][M+1];
    double Rho_bar[N+1][M+1];
    double Eta_bar[N+1][M+1];
    double X[M+1];
    double Y[N+1];
    double theta[N+1][M+1];
    double d_theta[N+1][M+1];
};
//g(theta)
double g(double theta){
    if (theta>=1) {
        return 1;
    }
    else
        return 0;
}

//calculate pressure from theta
double theta_to_pressure(double theta){
    return p_c+g(theta)*beta*log(theta);
}

//initialization
void init_param(struct Tribo *node){
    double hs,delta;
    double p,P;
    double x,y;
    double z;
    for (int i=0; i<=M; i++) {
        node->X[i] = -1.0+i*DX;
    }
    for (int j=0; j<=N; j++) {
        node->Y[j] = -1.0+j*DY;
    }
    node->H0 = h0/(b*b/R);
    for (int j=0 ; j<=N; j++) {
        for (int i=0; i<=M; i++) {
            //theta
            node->theta[j][i] = 1.0;
            node->d_theta[j][i] = 0.0;
            //pressure
            p = theta_to_pressure(node->theta[j][i]);
            node->P[j][i] = p/P_h;
            //film thickness
            x = -b+i*dx;
            y = -b+j*dy;
            hs = R-sqrt(pow(R, 2)-pow(x, 2)-pow(y, 2));
            delta = 0.0;
            node->HS[j][i] = hs/(b*b/R);
            node->Delta[j][i] = delta/(b*b/R);
            node->H[j][i] = (h0+hs+delta)/(b*b/R);
            //density
            P = node->P[j][i];
            node->Rho_bar[j][i] = 1+(0.6e-9*P_h*P)/(1+1.7e-9*P_h*P);
            //viscosity
            z = alpha_0/(5.1e-9*(log(eta_0)+9.67));
            node->Eta_bar[j][i] = exp((log(eta_0)+9.67) * (-1 + pow(1+P_h*P/1.98e8, z)));
        }
    }
    
}

//cal_epsilon = Rho_bar * H^3 / Eta_bar
double cal_epsilon(double Rho_bar, double H, double Eta_bar){
    return Rho_bar*pow(H, 3) / Eta_bar;
}
//g(theta-1)
double gTheta_m1(double theta){
    return g(theta)*(theta-1);
}

//calculate F(i,j)
double cal_F_ij(int i, int j, double Eta_bar[N+1][M+1], double Rho_bar[N+1][M+1], double H[N+1][M+1], double theta[N+1][M+1]){
    double A_ij,B_ij,C_ij;
    double epsilon_ij,epsilon_ip1j,epsilon_im1j,epsilon_ijm1,epsilon_ijp1;
    double Rho_H_theta_ip1j,Rho_H_theta_ij,Rho_H_theta_im1j;
    
    epsilon_ij   = cal_epsilon(Rho_bar[j][i],H[j][i],Eta_bar[j][i]);
    epsilon_ip1j = cal_epsilon(Rho_bar[j][i+1],H[j][i+1],Eta_bar[j][i+1]);
    epsilon_im1j = cal_epsilon(Rho_bar[j][i-1],H[j][i-1],Eta_bar[j][i-1]);
    epsilon_ijp1 = cal_epsilon(Rho_bar[j+1][i], H[j+1][i], Eta_bar[j+1][i]);
    epsilon_ijm1 = cal_epsilon(Rho_bar[j-1][i], H[j-1][i], Eta_bar[j-1][i]);
    
    Rho_H_theta_ij = Rho_bar[j][i]*H[j][i]*theta[j][i];
    Rho_H_theta_ip1j = Rho_bar[j][i+1]*H[j][i+1]*theta[j][i+1];
    Rho_H_theta_im1j = Rho_bar[j][i-1]*H[j][i-1]*theta[j][i-1];
    
    A_ij = 1/(2*DX) * ( (epsilon_ip1j+epsilon_ij)*gTheta_m1(theta[j][i+1]) - (epsilon_ip1j+2*epsilon_ij+epsilon_im1j)*gTheta_m1(theta[j][i]) + (epsilon_im1j+epsilon_ij)*gTheta_m1(theta[j][i-1]) );
    B_ij = 1/(2*DY) * ( (epsilon_ijp1+epsilon_ij)*gTheta_m1(theta[j+1][i]) - (epsilon_ijp1+2*epsilon_ij+epsilon_ijm1)*gTheta_m1(theta[j][i]) + (epsilon_ijm1+epsilon_ij)*gTheta_m1(theta[j-1][i]) );
    C_ij = ((1-beta_x)*(Rho_H_theta_ip1j-Rho_H_theta_ij) + (beta_x)*(Rho_H_theta_ij-Rho_H_theta_im1j)) /DX;
    
    return A_ij + k*k*B_ij - Psi*C_ij;
}

//calculate Ro,Et,Dm
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
//Th_ij_kl = d(theta_ij)/d(theta_kl)
double cal_Th_ij_kl(int i, int j, int k, int l, double theta[N+1][M+1]){
    if (i == k && j == l) return 1.0;
    return 0.0;
}
double cal_M_ij_kl(int i, int j, int k, int l, double Eta_bar[N+1][M+1], double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double H3_Eta_ij, Rho_H3_Eta2_ij, Rho_H2_Eta_ij;
    double Et_ij_kl,Rh_ij_kl, D_mn;
    int m,n;
    double dP_dtheta_kl;
    
    m = abs(k-i);
    n = abs(l-j);
    Rh_ij_kl = cal_Ro_ij_kl(i, j, k, l, P[j][i]);
    Et_ij_kl = cal_Et_ij_kl(i, j, k, l, P[j][i], Eta_bar[j][i]);
    D_mn = cal_D_mn(m, n);
    
    H3_Eta_ij = pow(H[j][i], 3)/Eta_bar[j][i];
    Rho_H2_Eta_ij = Rho_bar[j][i]*pow(H[j][i], 2)/Eta_bar[j][i];
    Rho_H3_Eta2_ij = Rho_bar[j][i]*pow(H[j][i], 3)/pow(Eta_bar[j][i], 2);
    
    dP_dtheta_kl = g(theta[l][k])*beta/( theta[l][k]* P_h);
    return dP_dtheta_kl * (H3_Eta_ij*Rh_ij_kl - Et_ij_kl*Rho_H3_Eta2_ij + 3*Rho_H2_Eta_ij * D_mn);
}
double cal_Nx_ij_kl(int i, int j, int k, int l, double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double Ro_ij_kl, D_mn,Th_ij_kl;
    double H_theta_ij, Rho_theta_ij;
    int m,n;
    double dP_dtheta_kl;
    
    m = abs(k-i);
    n = abs(l-j);
    H_theta_ij = H[j][i]*theta[j][i];
    Rho_theta_ij = Rho_bar[j][i]*theta[j][i];
    
    D_mn = cal_D_mn(m, n);
    Ro_ij_kl = cal_Ro_ij_kl(i,j, k, l, P[j][i]);
    Th_ij_kl = cal_Th_ij_kl(i, j, k, l, theta);
    
    dP_dtheta_kl = g(theta[l][k])*beta/( theta[l][k]* P_h);
    
    return (Rho_bar[j][i]*H[j][i])*Th_ij_kl + (Rho_bar[j][i]*theta[j][i])*D_mn*dP_dtheta_kl + (H[j][i]*theta[j][i])*Ro_ij_kl*dP_dtheta_kl;
}


//J(ij,ij) = dF(i,j)/dtheta(i,j)
double cal_J_ij_ij(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double epsilon_im1j,epsilon_ij,epsilon_ip1j,epsilon_ijm1,epsilon_ijp1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j, i, j, Eta_bar, Rho_bar, H, P, theta);
    M2 = cal_M_ij_kl(i+1, j, i, j, Eta_bar, Rho_bar, H, P, theta);
    M3 = cal_M_ij_kl(i  , j, i, j, Eta_bar, Rho_bar, H, P, theta);
    M4 = cal_M_ij_kl(i, j-1, i, j, Eta_bar, Rho_bar, H, P, theta);
    M5 = cal_M_ij_kl(i, j+1, i, j, Eta_bar, Rho_bar, H, P, theta);
    
    Nx1 = cal_Nx_ij_kl(i-1, j, i, j, Rho_bar, H, P, theta);
    Nx2 = cal_Nx_ij_kl(i  , j, i, j, Rho_bar, H, P, theta);
    Nx3 = cal_Nx_ij_kl(i+1, j, i, j, Rho_bar, H, P, theta);
    
    epsilon_ij   = cal_epsilon(Rho_bar[j][i]  , H[j][i]  , Eta_bar[j][i]  );
    epsilon_im1j = cal_epsilon(Rho_bar[j][i-1], H[j][i-1], Eta_bar[j][i-1]);
    epsilon_ip1j = cal_epsilon(Rho_bar[j][i+1], H[j][i+1], Eta_bar[j][i+1]);
    epsilon_ijm1 = cal_epsilon(Rho_bar[j-1][i], H[j-1][i], Eta_bar[j-1][i]);
    epsilon_ijp1 = cal_epsilon(Rho_bar[j+1][i], H[j+1][i], Eta_bar[j+1][i]);
    
    A1 = 1  /(2*DX*DX) * (gTheta_m1(theta[j][i-1])*(M3+M1) - gTheta_m1(theta[j][i])*(M2+M3*2+M1) + gTheta_m1(theta[j][i+1])*(M2+M3) - g(theta[j][i])*(epsilon_im1j+epsilon_ip1j+2*epsilon_ij));
    A2 = k*k/(2*DY*DY) * (gTheta_m1(theta[j-1][i])*(M4+M3) - gTheta_m1(theta[j][i])*(M4+M3*2+M5) + gTheta_m1(theta[j+1][i])*(M5+M3) - g(theta[j][i])*(epsilon_ijm1+epsilon_ijp1+2*epsilon_ij));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C = 0.0;
    
    return A1+A2-B1-B2-C;
    
}

//J(ij,i+1j) = dF(i,j)/dtheta(i+1,j)
double cal_J_ij_ip1j(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double epsilon_im1j,epsilon_ij,epsilon_ip1j,epsilon_ijp1,epsilon_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j,   i+1, j, Eta_bar, Rho_bar, H, P, theta);
    M2 = cal_M_ij_kl(i+1, j,   i+1, j, Eta_bar, Rho_bar, H, P, theta);
    M3 = cal_M_ij_kl(i,   j,   i+1, j, Eta_bar, Rho_bar, H, P, theta);
    M4 = cal_M_ij_kl(i,   j-1, i+1, j, Eta_bar, Rho_bar, H, P, theta);
    M5 = cal_M_ij_kl(i,   j+1, i+1, j, Eta_bar, Rho_bar, H, P, theta);
    
    Nx1 = cal_Nx_ij_kl(i-1, j, i+1, j, Rho_bar, H, P, theta);
    Nx2 = cal_Nx_ij_kl(i,   j, i+1, j, Rho_bar, H, P, theta);
    Nx3 = cal_Nx_ij_kl(i+1, j, i+1, j, Rho_bar, H, P, theta);
    
    epsilon_ij   = cal_epsilon(Rho_bar[j][i],   H[j][i]  , Eta_bar[j][i]);
    epsilon_im1j = cal_epsilon(Rho_bar[j][i-1], H[j][i-1], Eta_bar[j][i-1]);
    epsilon_ip1j = cal_epsilon(Rho_bar[j][i+1], H[j][i+1], Eta_bar[j][i+1]);
    epsilon_ijm1 = cal_epsilon(Rho_bar[j-1][i], H[j-1][i], Eta_bar[j-1][i]);
    epsilon_ijp1 = cal_epsilon(Rho_bar[j+1][i], H[j+1][i], Eta_bar[j+1][i]);
    
    A1 = 1  /(2*DX*DX) * (gTheta_m1(theta[j][i-1])*(M3+M1) - gTheta_m1(theta[j][i])*(M2+M3*2+M1) + gTheta_m1(theta[j][i+1])*(M2+M3) + g(theta[j][i+1])*(epsilon_ip1j+epsilon_ij));
    A2 = k*k/(2*DY*DY) * (gTheta_m1(theta[j-1][i])*(M4+M3) - gTheta_m1(theta[j][i])*(M4+M3*2+M5) + gTheta_m1(theta[j+1][i])*(M5+M3));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C = 0.0;
    
    return A1+A2-B1-B2-C;
}

//J(ij,i-1j) = dF(i,j)/dtheta(i-1,j)
double cal_J_ij_im1j(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double epsilon_im1j,epsilon_ij,epsilon_ip1j,epsilon_ijp1,epsilon_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j,   i-1, j, Eta_bar, Rho_bar, H, P, theta);
    M2 = cal_M_ij_kl(i+1, j,   i-1, j, Eta_bar, Rho_bar, H, P, theta);
    M3 = cal_M_ij_kl(i,   j,   i-1, j, Eta_bar, Rho_bar, H, P, theta);
    M4 = cal_M_ij_kl(i,   j-1, i-1, j, Eta_bar, Rho_bar, H, P, theta);
    M5 = cal_M_ij_kl(i,   j+1, i-1, j, Eta_bar, Rho_bar, H, P, theta);
    
    Nx1 = cal_Nx_ij_kl(i-1, j, i-1, j, Rho_bar, H, P, theta);
    Nx2 = cal_Nx_ij_kl(i,   j, i-1, j, Rho_bar, H, P, theta);
    Nx3 = cal_Nx_ij_kl(i+1, j, i-1, j, Rho_bar, H, P, theta);
    
    epsilon_ij   = cal_epsilon(Rho_bar[j][i],   H[j][i],   Eta_bar[j][i]  );
    epsilon_im1j = cal_epsilon(Rho_bar[j][i-1], H[j][i-1], Eta_bar[j][i-1]);
    epsilon_ip1j = cal_epsilon(Rho_bar[j][i+1], H[j][i+1], Eta_bar[j][i+1]);
    epsilon_ijm1 = cal_epsilon(Rho_bar[j-1][i], H[j-1][i], Eta_bar[j-1][i]);
    epsilon_ijp1 = cal_epsilon(Rho_bar[j+1][i], H[j+1][i], Eta_bar[j+1][i]);
    
    A1 = 1  /(2*DX*DX) * (gTheta_m1(theta[j][i-1])*(M3+M1) - gTheta_m1(theta[j][i])*(M2+M3*2+M1) + gTheta_m1(theta[j][i+1])*(M2+M3) + g(theta[j][i-1])*(epsilon_im1j+epsilon_ij));
    A2 = k*k/(2*DY*DY) * (gTheta_m1(theta[j-1][i])*(M4+M3) - gTheta_m1(theta[j][i])*(M4+M3*2+M5) + gTheta_m1(theta[j+1][i])*(M5+M3));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C = 0.0;
    
    return A1+A2-B1-B2-C;
}
//J(ij,ij+1) = dF(i,j)/dtheta(i,j+1)
double cal_J_ij_ijp1(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double epsilon_im1j,epsilon_ij,epsilon_ip1j,epsilon_ijp1,epsilon_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j,   i, j+1, Eta_bar, Rho_bar, H, P, theta);
    M2 = cal_M_ij_kl(i+1, j,   i, j+1, Eta_bar, Rho_bar, H, P, theta);
    M3 = cal_M_ij_kl(i,   j,   i, j+1, Eta_bar, Rho_bar, H, P, theta);
    M4 = cal_M_ij_kl(i,   j-1, i, j+1, Eta_bar, Rho_bar, H, P, theta);
    M5 = cal_M_ij_kl(i,   j+1, i, j+1, Eta_bar, Rho_bar, H, P, theta);
    
    Nx1 = cal_Nx_ij_kl(i-1, j, i, j+1, Rho_bar, H, P, theta);
    Nx2 = cal_Nx_ij_kl(i,   j, i, j+1, Rho_bar, H, P, theta);
    Nx3 = cal_Nx_ij_kl(i+1, j, i, j+1, Rho_bar, H, P, theta);
    
    epsilon_ij   = cal_epsilon(Rho_bar[j][i],   H[j][i],   Eta_bar[j][i]  );
    epsilon_im1j = cal_epsilon(Rho_bar[j][i-1], H[j][i-1], Eta_bar[j][i-1]);
    epsilon_ip1j = cal_epsilon(Rho_bar[j][i+1], H[j][i+1], Eta_bar[j][i+1]);
    epsilon_ijm1 = cal_epsilon(Rho_bar[j-1][i], H[j-1][i], Eta_bar[j-1][i]);
    epsilon_ijp1 = cal_epsilon(Rho_bar[j+1][i], H[j+1][i], Eta_bar[j+1][i]);
    
    A1 = 1  /(2*DX*DX) * (gTheta_m1(theta[j][i-1])*(M3+M1) - gTheta_m1(theta[j][i])*(M2+M3*2+M1) + gTheta_m1(theta[j][i+1])*(M2+M3));
    A2 = k*k/(2*DY*DY) * (gTheta_m1(theta[j-1][i])*(M4+M3) - gTheta_m1(theta[j][i])*(M4+M3*2+M5) + gTheta_m1(theta[j+1][i])*(M5+M3) + g(theta[j+1][i])*(epsilon_ip1j+epsilon_ij));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C = 0.0;
    
    return A1+A2-B1-B2-C;
}
//J(ij,ij-1) = dF(i,j)/dtheta(i,j-1)
double cal_J_ij_ijm1(int i,int j,double Eta_bar[N+1][M+1],double Rho_bar[N+1][M+1], double H[N+1][M+1], double P[N+1][M+1], double theta[N+1][M+1]){
    double M1,M2,M3,M4,M5;
    double Nx1,Nx2,Nx3;
    double epsilon_im1j,epsilon_ij,epsilon_ip1j,epsilon_ijp1,epsilon_ijm1;
    double A1,A2,B1,B2,C;
    
    M1 = cal_M_ij_kl(i-1, j,   i, j-1, Eta_bar, Rho_bar, H, P, theta);
    M2 = cal_M_ij_kl(i+1, j,   i, j-1, Eta_bar, Rho_bar, H, P, theta);
    M3 = cal_M_ij_kl(i,   j,   i, j-1, Eta_bar, Rho_bar, H, P, theta);
    M4 = cal_M_ij_kl(i,   j-1, i, j-1, Eta_bar, Rho_bar, H, P, theta);
    M5 = cal_M_ij_kl(i,   j+1, i, j-1, Eta_bar, Rho_bar, H, P, theta);
    
    Nx1 = cal_Nx_ij_kl(i-1, j, i, j-1, Rho_bar, H, P, theta);
    Nx2 = cal_Nx_ij_kl(i,   j, i, j-1, Rho_bar, H, P, theta);
    Nx3 = cal_Nx_ij_kl(i+1, j, i, j-1, Rho_bar, H, P, theta);
    
    epsilon_ij   = cal_epsilon(Rho_bar[j][i],   H[j][i],   Eta_bar[j][i]  );
    epsilon_im1j = cal_epsilon(Rho_bar[j][i-1], H[j][i-1], Eta_bar[j][i-1]);
    epsilon_ip1j = cal_epsilon(Rho_bar[j][i+1], H[j][i+1], Eta_bar[j][i+1]);
    epsilon_ijm1 = cal_epsilon(Rho_bar[j-1][i], H[j-1][i], Eta_bar[j-1][i]);
    epsilon_ijp1 = cal_epsilon(Rho_bar[j+1][i], H[j+1][i], Eta_bar[j+1][i]);
    
    A1 = 1  /(2*DX*DX) * (gTheta_m1(theta[j][i-1])*(M3+M1) - gTheta_m1(theta[j][i])*(M2+M3*2+M1) + gTheta_m1(theta[j][i+1])*(M2+M3));
    A2 = k*k/(2*DY*DY) * (gTheta_m1(theta[j-1][i])*(M4+M3) - gTheta_m1(theta[j][i])*(M4+M3*2+M5) + gTheta_m1(theta[j+1][i])*(M5+M3) + g(theta[j-1][i])*(epsilon_im1j+epsilon_ij));
    B1 = Psi/DX*((1-beta_x)*(Nx3-Nx2) + beta_x*(Nx2-Nx1));
    B2 = 0.0;
    C = 0.0;
    
    return A1+A2-B1-B2-C;
}

void update_param(struct Tribo *node){
    double F,J1,J2,J3,J4,J5;
    double delta_kl;
    double z;
    E_dash = 1/(1/E1+1/E2);
    z = alpha_0/(5.1e-9*(log(eta_0)+9.67));
    for (int j=1; j < N ; j++) {
        printf("j=%d\n",j);
        for (int i = 1; i<M; i++) {
            //calculate dtheta
            F  = cal_F_ij(i, j, node->Eta_bar, node->Rho_bar, node->H, node->theta); //F(i,j)
            J1 = cal_J_ij_im1j(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P,node->theta); //J(i-1,j)
            J2 = cal_J_ij_ip1j(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P,node->theta); //J(i+1,j)
            J3 = cal_J_ij_ij(  i, j, node->Eta_bar, node->Rho_bar, node->H, node->P,node->theta); //J(i,j)
            J4 = cal_J_ij_ijm1(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P,node->theta); //J(i,j-1)
            J5 = cal_J_ij_ijp1(i, j, node->Eta_bar, node->Rho_bar, node->H, node->P,node->theta); //J(i,j+1)
            node->d_theta[j][i] = -(F + J1*node->d_theta[j][i-1] + J2*node->d_theta[j][i+1] + J4*node->d_theta[j-1][i] + J5*node->d_theta[j+1][i]) / J3;
            printf("%lf\t",node->d_theta[j][i]);
            //update theta
            node->theta[j][i] = node->theta[j][i] + omega*node->d_theta[j][i];
            if (node->theta[j][i] <= 0.0) {
                node->theta[j][i] = 1.0e-6;
            }
            //update pressure
            node->P[j][i]= (g(node->theta[j][i])*beta*log(node->theta[j][i])+p_c)/P_h;
            //update viscosity
            node->Eta_bar[j][i] = exp((log(eta_0)+9.67) * (-1 + pow(1+P_h*node->P[j][i]/1.98e8, z)));
            //update density
            node->Rho_bar[j][i] = 1+(0.6e-9*P_h*node->P[j][i])/(1+1.7e-9*P_h*node->P[j][i]);
            //update film thickness
            for (int l= 0; l <= N; l++) {
                for (int k = 0; k <= M; k++ ) {
                    delta_kl = cal_delta_kl(k, l, node->P);
                    node->Delta[l][k] = delta_kl/(b*b/R);
                    node->H[l][k] = node->H0+node->HS[l][k]+node->Delta[l][k];
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
//print the parameter
void printf_param(double matrix[N+1][M+1]){
    for (int j=0; j<=N; j++) {
        for (int i=0; i<=M; i++) {
            printf("%.5lf\t",matrix[j][i]);
        }printf("\n");
    }
}

int main(int argc, const char * argv[]) {
    struct Tribo node;
    double theta_pervious[N+1][M+1], theta_new[N+1][M+1];
    int loop = 0;
    
    dx = 2*b/M;
    dy = 2*a/N;
    DX = dx / b;
    DY = dy / a;
    
    beta_bar = beta*R/(eta_0*u_avg);
    Psi = 12*pow(R,3)/(pow(b,3)*beta_bar);
    //initialize the nodes and all the parameter
    init_param(&node);

    double diff_theta, sum_P = 0.0;
    double pMax = 0.0, hMin = 1e10;
    //do-loop
    do{
        ++loop;
        diff_theta = 0.0;
        printf("======================loop: %d-th\t======================\n",loop);
        //save [k-1] th data
        copy_array(theta_pervious, node.theta);
        update_param(&node);
        //save [k] th data
        copy_array(theta_new, node.theta);
        for (int j=0; j<=N; j++) {
            for (int i=0; i<=M; i++) {
                diff_theta+=sqrt(fabs(theta_new[j][i]-theta_pervious[j][i]));
            }
        }
        diff_theta =diff_theta/M/N;
        printf("the error is :\t%lf\n",diff_theta);
        printf_param(node.theta);
    }while(diff_theta > Err_thetaMin || loop < 100);
    return 0;
}

