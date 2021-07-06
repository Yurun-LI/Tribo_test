#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI acos(-1)
#define R 7.5e-3
#define Omega 1.0e-3
#define M 100             //maximum node in x-axis
#define N 100             //maximum node in y-axis
#define Err_thetaMin 1.0e-7 //convergence criteria

double uA = 0.0209, uB = 0.0; // velocity 100rpm [m/s] in x-
double u_avg = 0.0209 / 2;    // average velocity [m/s]
double v = 0.0;               // velocity in y-axis
double a = 4e-6 / 2.0;        //[m]
double b = 4e-6 / 2.0;        //[m]
double k = 1.0;               // b/a
double U = 1.0, V = 0.0;      // non-dimensional velocity

double rho_0 = 1.82e3;                                           // [kg/m3]
double eta_0 = 0.055;                                            //[Pa.s]
double alpha_0 = 50.0e-9;                                        // [Pa-1]
double p_h = 1.01e6;                                             // [Pa_s]
double p_c = 0.02e6;                                             // [Pa_s]
double p_0 = 1.98e8;                                             // [Pa_s]
double E1 = 150.0e9, E2 = 110.0e9, E_dash, v1 = 0.25, v2 = 0.28; //Young's molecule [Pa]
double beta = 1.7e9;                                             // bulk modules of PFPE

double h_init = 1.7e-9; //initial film thickness
double h0 = 4e-9;       //clearance [m]

double dx, dy, DX, DY;
double beta_x = 0.5, beta_y = 0.5;
double Psi, beta_bar;

FILE *fp1, *fp2, *fp3;

struct Nodes
{
    //theta
    double theta[N + 1][M + 1];
    double d_theta[N + 1][M + 1];
    //pressure
    double P[N + 1][M + 1];
    //thickness
    double H[N + 1][M + 1];
    double HS[N + 1][M + 1];
    double H0;
    double Delta[N + 1][M + 1];
    //density and viscosity
    double Rho_c_bar[N + 1][M + 1];
    double Eta_bar[N + 1][M + 1];
    //locations
    double X[M + 1];
    double Y[N + 1];
};

//g(theta)
double g(double theta)
{
    if (theta >= 1)
        return 1.0;
    return 0.0;
}

//p(theta)
double pressure(double theta)
{
    return p_c + g(theta) * beta * log(theta);
}

//rho(pressure)
double rho_ij(double p)
{
    return rho_0 * (1 + (0.6 * 10e-9 * p) / (1 + 1.7 * 10e-9 * p));
}

//eta(pressure) --Roelands equation
double eta_ij(double p)
{
    double Z, alpha;
    Z = alpha_0 / (5.1e-9 * (log(eta_0) + 9.67));
    alpha = (log(eta_0) + 9.67) * (pow((1 + p / 1.98e8), Z) - 1);
    return eta_0 * exp(alpha);
}

//calculate Rho_c_bar * H^3 / Eta_bar  -> Epsilon_bar [i,j]
double cal_Epsilon_bar(double Rho_c_bar_ij, double H_ij, double Eta_bar_ij)
{
    return Rho_c_bar_ij * pow(H_ij, 3) / Eta_bar_ij;
}

//calculate g*(theta-1) [i,j]
double g_theta_m1(double theta_ij)
{
    return g(theta_ij) * (theta_ij - 1);
}

//calculate F(i,j)
double cal_F_ij(int i, int j, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double Eta_bar[M + 1][N + 1], double theta[N + 1][M + 1])
{
    double Ax_ij, Ay_ij, Bx_ij;
    double Epsilon_ij, Epsilon_ip1j, Epsilon_im1j, Epsilon_ijm1, Epsilon_ijp1;
    double Rho_H_theta_ip1j, Rho_H_theta_ij, Rho_H_theta_im1j, Rho_H_theta_ijp1, Rho_H_theta_ijm1;

    //Epsilon
    Epsilon_im1j = cal_Epsilon_bar(Rho_c_bar[j][i - 1], H[j][i - 1], Eta_bar[j][i - 1]);
    Epsilon_ip1j = cal_Epsilon_bar(Rho_c_bar[j][i + 1], H[j][i + 1], Eta_bar[j][i + 1]);
    Epsilon_ij = cal_Epsilon_bar(Rho_c_bar[j][i], H[j][i], Eta_bar[j][i]);
    Epsilon_ijm1 = cal_Epsilon_bar(Rho_c_bar[j - 1][i], H[j - 1][i], Eta_bar[j - 1][i]);
    Epsilon_ijp1 = cal_Epsilon_bar(Rho_c_bar[j + 1][i], H[j + 1][i], Eta_bar[j + 1][i]);

    //Rho_H_theta
    Rho_H_theta_im1j = Rho_c_bar[j][i - 1] * H[j][i - 1] * theta[j][i - 1];
    Rho_H_theta_ip1j = Rho_c_bar[j][i + 1] * H[j][i + 1] * theta[j][i + 1];
    Rho_H_theta_ij = Rho_c_bar[j][i] * H[j][i] * theta[j][i];

    Ax_ij = 1 / (2 * DX * DX) * (g_theta_m1(theta[j][i + 1]) * (Epsilon_ip1j + Epsilon_ij) - g_theta_m1(theta[j][i]) * (Epsilon_ip1j + Epsilon_ij * 2 + Epsilon_im1j) + g_theta_m1(theta[j][i - 1]) * (Epsilon_im1j + Epsilon_ij));
    Ay_ij = 1 / (2 * DY * DY) * (g_theta_m1(theta[j + 1][i]) * (Epsilon_ijp1 + Epsilon_ij) - g_theta_m1(theta[j][i]) * (Epsilon_ijp1 + Epsilon_ij * 2 + Epsilon_ijm1) + g_theta_m1(theta[j - 1][i]) * (Epsilon_ijm1 + Epsilon_ij));
    Bx_ij = (1 / DX) * ((1 - beta_x) * (Rho_H_theta_ip1j - Rho_H_theta_ij) + (beta_x) * (Rho_H_theta_ij - Rho_H_theta_im1j));

    return Ax_ij + k * k * Ay_ij - Psi * Bx_ij;
}
double dtheta_dgtheta(int i, int j, int k, int l, double theta[N + 1][M + 1])
{
    if ((theta[j][i] >= 1) && (i == k) && (j == l))
        return 1.0;
    return 0.0;
}
double cal_N_ij_kl(int i, int j, int k, int l, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double theta[N + 1][M + 1])
{
    return Rho_c_bar[j][i] * H[j][i] * dtheta_dgtheta(i, j, k, l, theta);
}
//Jacobian, J terms
//J0 = dF[i,j]/d(g_theta)[i+1,j]
double cal_J0(int i, int j, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double Eta_bar[N + 1][M + 1], double theta[N + 1][M + 1])
{
    double Ax0, Ay0, B0;
    double M1, M2, M3, M4, M5;
    double N1, N2, N3;
    double Epsilon_ij, Epsilon_ip1j, Epsilon_im1j, Epsilon_ijm1, Epsilon_ijp1;

    Epsilon_im1j = cal_Epsilon_bar(Rho_c_bar[j][i - 1], H[j][i - 1], Eta_bar[j][i - 1]);
    Epsilon_ip1j = cal_Epsilon_bar(Rho_c_bar[j][i + 1], H[j][i + 1], Eta_bar[j][i + 1]);
    Epsilon_ij = cal_Epsilon_bar(Rho_c_bar[j][i], H[j][i], Eta_bar[j][i]);
    Epsilon_ijm1 = cal_Epsilon_bar(Rho_c_bar[j - 1][i], H[j - 1][i], Eta_bar[j - 1][i]);
    Epsilon_ijp1 = cal_Epsilon_bar(Rho_c_bar[j + 1][i], H[j + 1][i], Eta_bar[j + 1][i]);

    N1 = cal_N_ij_kl(i - 1, j, i + 1, j, Rho_c_bar, H, theta);
    N2 = cal_N_ij_kl(i, j, i + 1, j, Rho_c_bar, H, theta);
    N3 = cal_N_ij_kl(i + 1, j, i + 1, j, Rho_c_bar, H, theta);

    Ax0 = 1 / (2 * DX * DX) * ((Epsilon_ip1j + Epsilon_ij));
    Ay0 = 0;
    B0 = (1 / DX) * ((1 - beta_x) * (N3 - N2) + (beta_x) * (N2 - N1));

    return Ax0 + k * k * Ay0 - Psi * B0;
}

//J1 = dF[i,j]/d(g_theta)[i-1,j]
double cal_J1(int i, int j, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double Eta_bar[N + 1][M + 1], double theta[N + 1][M + 1])
{
    double Ax0, Ay0, B0;
    double M1, M2, M3, M4, M5;
    double N1, N2, N3;
    double Epsilon_ij, Epsilon_ip1j, Epsilon_im1j, Epsilon_ijm1, Epsilon_ijp1;

    Epsilon_im1j = cal_Epsilon_bar(Rho_c_bar[j][i - 1], H[j][i - 1], Eta_bar[j][i - 1]);
    Epsilon_ip1j = cal_Epsilon_bar(Rho_c_bar[j][i + 1], H[j][i + 1], Eta_bar[j][i + 1]);
    Epsilon_ij = cal_Epsilon_bar(Rho_c_bar[j][i], H[j][i], Eta_bar[j][i]);
    Epsilon_ijm1 = cal_Epsilon_bar(Rho_c_bar[j - 1][i], H[j - 1][i], Eta_bar[j - 1][i]);
    Epsilon_ijp1 = cal_Epsilon_bar(Rho_c_bar[j + 1][i], H[j + 1][i], Eta_bar[j + 1][i]);

    N1 = cal_N_ij_kl(i - 1, j, i - 1, j, Rho_c_bar, H, theta);
    N2 = cal_N_ij_kl(i, j, i - 1, j, Rho_c_bar, H, theta);
    N3 = cal_N_ij_kl(i + 1, j, i - 1, j, Rho_c_bar, H, theta);

    Ax0 = 1 / (2 * DX * DX) * ((Epsilon_im1j + Epsilon_ij));
    Ay0 = 0;
    B0 = (1 / DX) * ((1 - beta_x) * (N3 - N2) + (beta_x) * (N2 - N1));

    return Ax0 + k * k * Ay0 - Psi * B0;
}

//J2 = dF[i,j]/d(g_theta)[i,j+1]
double cal_J2(int i, int j, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double Eta_bar[N + 1][M + 1], double theta[N + 1][M + 1])
{
    double Ax0, Ay0, B0;
    double M1, M2, M3, M4, M5;
    double N1, N2, N3;
    double Epsilon_ij, Epsilon_ip1j, Epsilon_im1j, Epsilon_ijm1, Epsilon_ijp1;

    Epsilon_im1j = cal_Epsilon_bar(Rho_c_bar[j][i - 1], H[j][i - 1], Eta_bar[j][i - 1]);
    Epsilon_ip1j = cal_Epsilon_bar(Rho_c_bar[j][i + 1], H[j][i + 1], Eta_bar[j][i + 1]);
    Epsilon_ij = cal_Epsilon_bar(Rho_c_bar[j][i], H[j][i], Eta_bar[j][i]);
    Epsilon_ijm1 = cal_Epsilon_bar(Rho_c_bar[j - 1][i], H[j - 1][i], Eta_bar[j - 1][i]);
    Epsilon_ijp1 = cal_Epsilon_bar(Rho_c_bar[j + 1][i], H[j + 1][i], Eta_bar[j + 1][i]);

    N1 = cal_N_ij_kl(i - 1, j, i, j + 1, Rho_c_bar, H, theta);
    N2 = cal_N_ij_kl(i, j, i, j + 1, Rho_c_bar, H, theta);
    N3 = cal_N_ij_kl(i + 1, j, i, j + 1, Rho_c_bar, H, theta);

    Ax0 = 0.0;
    Ay0 = 1 / (2 * DY * DY) * ((Epsilon_ijp1 + Epsilon_ij));
    B0 = (1 / DX) * ((1 - beta_x) * (N3 - N2) + (beta_x) * (N2 - N1));

    return Ax0 + k * k * Ay0 - Psi * B0;
}

//J3 = dF[i,j]/d(g_theta)[i,j-1]
double cal_J3(int i, int j, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double Eta_bar[N + 1][M + 1], double theta[N + 1][M + 1])
{
    double Ax0, Ay0, B0;
    double M1, M2, M3, M4, M5;
    double N1, N2, N3;
    double Epsilon_ij, Epsilon_ip1j, Epsilon_im1j, Epsilon_ijm1, Epsilon_ijp1;

    Epsilon_im1j = cal_Epsilon_bar(Rho_c_bar[j][i - 1], H[j][i - 1], Eta_bar[j][i - 1]);
    Epsilon_ip1j = cal_Epsilon_bar(Rho_c_bar[j][i + 1], H[j][i + 1], Eta_bar[j][i + 1]);
    Epsilon_ij = cal_Epsilon_bar(Rho_c_bar[j][i], H[j][i], Eta_bar[j][i]);
    Epsilon_ijm1 = cal_Epsilon_bar(Rho_c_bar[j - 1][i], H[j - 1][i], Eta_bar[j - 1][i]);
    Epsilon_ijp1 = cal_Epsilon_bar(Rho_c_bar[j + 1][i], H[j + 1][i], Eta_bar[j + 1][i]);

    N1 = cal_N_ij_kl(i - 1, j, i, j - 1, Rho_c_bar, H, theta);
    N2 = cal_N_ij_kl(i, j, i, j - 1, Rho_c_bar, H, theta);
    N3 = cal_N_ij_kl(i + 1, j, i, j - 1, Rho_c_bar, H, theta);

    Ax0 = 0.0;
    Ay0 = 1 / (2 * DY * DY) * ((Epsilon_ijm1 + Epsilon_ij));
    B0 = (1 / DX) * ((1 - beta_x) * (N3 - N2) + (beta_x) * (N2 - N1));

    return Ax0 + k * k * Ay0 - Psi * B0;
}

//J4 = dF[i,j]/d(g_theta)[i,j]
double cal_J4(int i, int j, double Rho_c_bar[N + 1][M + 1], double H[N + 1][M + 1], double Eta_bar[N + 1][M + 1], double theta[N + 1][M + 1])
{
    double Ax0, Ay0, B0;
    double M1, M2, M3, M4, M5;
    double N1, N2, N3;
    double Epsilon_ij, Epsilon_ip1j, Epsilon_im1j, Epsilon_ijm1, Epsilon_ijp1;

    Epsilon_im1j = cal_Epsilon_bar(Rho_c_bar[j][i - 1], H[j][i - 1], Eta_bar[j][i - 1]);
    Epsilon_ip1j = cal_Epsilon_bar(Rho_c_bar[j][i + 1], H[j][i + 1], Eta_bar[j][i + 1]);
    Epsilon_ij = cal_Epsilon_bar(Rho_c_bar[j][i], H[j][i], Eta_bar[j][i]);
    Epsilon_ijm1 = cal_Epsilon_bar(Rho_c_bar[j - 1][i], H[j - 1][i], Eta_bar[j - 1][i]);
    Epsilon_ijp1 = cal_Epsilon_bar(Rho_c_bar[j + 1][i], H[j + 1][i], Eta_bar[j + 1][i]);

    N1 = cal_N_ij_kl(i - 1, j, i, j, Rho_c_bar, H, theta);
    N2 = cal_N_ij_kl(i, j, i, j, Rho_c_bar, H, theta);
    N3 = cal_N_ij_kl(i + 1, j, i, j, Rho_c_bar, H, theta);

    Ax0 = 1 / (2 * DX * DX) * (0.0 - (Epsilon_ip1j + 2 * Epsilon_ij + Epsilon_im1j));
    Ay0 = 1 / (2 * DY * DY) * (0.0 - (Epsilon_ijp1 + 2 * Epsilon_ij + Epsilon_ijm1));
    B0 = (1 / DX) * ((1 - beta_x) * (N3 - N2) + (beta_x) * (N2 - N1));

    return Ax0 + k * k * Ay0 - Psi * B0;
}

//initialization
void init_param(struct Nodes *node)
{
    double hs, delta;
    double p;
    double rho, rho_c;
    double eta;
    double x, y;

    for (int i = 0; i <= M; i++)
        node->X[i] = -1.0 + i * DX;
    for (int j = 0; j <= N; j++)
        node->Y[j] = -1.0 + j * DY;

    node->H0 = h0 / (b * b / R);

    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= M; i++)
        {
            //theta
            node->theta[j][i] = 1.0;
            node->d_theta[j][i] = 0.0;

            //theta->pressure->dimensionless pressure
            p = pressure(node->theta[j][i]);
            node->P[j][i] = p / p_h; //dimensionless pressure

            //pressure->density->dimensionless density
            rho = rho_ij(p);
            rho_c = rho / node->theta[j][i];
            node->Rho_c_bar[j][i] = rho_c / rho_0; //dimensionless density

            //pressure->viscosity->dimensionless viscosity
            eta = eta_ij(p);
            node->Eta_bar[j][i] = eta / eta_0; //dimensionless viscosity

            //pressure->film thickness -> dimensionless film thickness
            x = -b + i * dx;
            y = -b + j * dy;
            hs = R - sqrt(pow(R, 2) - pow(x, 2) - pow(y, 2));
            delta = 0.0;
            node->HS[j][i] = hs / (b * b / R);
            node->Delta[j][i] = delta / (b * b / R);
            node->H[j][i] = node->H0 + node->HS[j][i] + node->Delta[j][i];
        }
    }
    printf("all the nodes are initialized!\n");
}

//copy array
void copy_array(double new_array[N + 1][M + 1], double array[N + 1][M + 1])
{
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= M; i++)
        {
            new_array[j][i] = array[j][i];
        }
    }
}

//print the parameter
void printf_param(double matrix[N + 1][M + 1])
{
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= M; i++)
        {
            printf("%.5lf\t", matrix[j][i]);
        }
        printf("\n");
    }
}

//update parameters
void update_param(struct Nodes *node)
{
    double F, J0, J1, J2, J3, J4;
    double p;
    for (int j = 1; j < N; j++)
    {
        for (int i = 1; i < M; i++)
        {
            F = cal_F_ij(i, j, node->Rho_c_bar, node->H, node->Eta_bar, node->theta);
            J0 = cal_J0(i, j, node->Rho_c_bar, node->H, node->Eta_bar, node->theta);
            J1 = cal_J1(i, j, node->Rho_c_bar, node->H, node->Eta_bar, node->theta);
            J2 = cal_J2(i, j, node->Rho_c_bar, node->H, node->Eta_bar, node->theta);
            J3 = cal_J3(i, j, node->Rho_c_bar, node->H, node->Eta_bar, node->theta);
            J4 = cal_J4(i, j, node->Rho_c_bar, node->H, node->Eta_bar, node->theta);

            node->d_theta[j][i] = -(F + J0 * node->d_theta[j][i + 1] + J1 * node->d_theta[j][i - 1] + J2 * node->d_theta[j + 1][i] + J3 * node->d_theta[j - 1][i]) / J4;
            
            //update theta
            node->theta[j][i] = node->theta[j][i] + Omega * node->d_theta[j][i];
            if (node->theta[j][i] <= 0.0)
            {
                node->theta[j][i] = 1.0e-6;
            }
            //update pressure
            p = pressure(node->theta[j][i]);
            node->P[j][i] = p / p_h;
        }
    }
}

int main(int argc, char *argv[])
{
    struct Nodes node;
    double theta_pervious[N + 1][M + 1], theta_new[N + 1][M + 1];
    int theta_loop = 0;

    double error_theta = 0.0;

    dx = 2 * b / M;
    dy = 2 * a / N;
    DX = dx / b;
    DY = dy / a;

    beta_bar = beta * R / (eta_0 * u_avg);
    Psi = 12 * pow(R, 3) / (pow(b, 3) * beta_bar);

    //initialize nodes
    init_param(&node);

    do //Gauss-Seidal iteration
    {
        ++theta_loop;
        double diff_theta = 0.0, sum_theta = 0.0;
        printf("---------------------\tloop: %d-th\t---------------------\n", theta_loop);
        //save the (k-1)th data
        copy_array(theta_pervious, node.theta);
        update_param(&node);
        //save the k-th data
        copy_array(theta_new, node.theta);
        for (int j = 0; j <= N; j++)
        {
            for (int i = 0; i <= M; i++)
            {
                diff_theta += fabs(theta_new[j][i] - theta_pervious[j][i]);
                sum_theta += node.theta[j][i];
            }
        }
        error_theta = diff_theta / sum_theta;
//        printf_param(node.theta);
        printf("the error is :\t%.10lf\n\n", error_theta);
    } while (error_theta > Err_thetaMin || theta_loop <= 1000);
    printf_param(node.theta);
    printf("\n");
    printf_param(node.P);
    //save data
    fp1 = fopen("/Users/liyurun/Desktop/tribometer_data/data_for_20210706/node_100_constant_theta.csv", "w+");
    for (int j = 0; j<N+1; j++){
        for (int i = 0; i< M+1; i++){
            fprintf(fp1, "%.5lf%s",node.theta[j][i], //output pressure [Pa.s]
                    (i<M+1-1?",":""));
        }
        fprintf(fp1,"\n");
    }
    fclose(fp1);

    fp2 = fopen("/Users/liyurun/Desktop/tribometer_data/data_for_20210706/node_100_constant_pressure.csv", "w+");
    for (int j = 0; j<N+1; j++){
        for (int i = 0; i< M+1; i++){
            fprintf(fp2, "%.10lf%s",node.P[j][i]*p_h, //output film thickness [m]
                    (i<M+1-1?",":""));
        }
        fprintf(fp2,"\n");
    }
    fclose(fp2);
    
    return 0;
}

