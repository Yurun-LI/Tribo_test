#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define pi acos(-1) //pi in radian

#define MAX_x 100 //max range in x-direction (point)

#define boundary 2.0 //boundary condition

#define beta 1.0 //beta in finite difference method

#define omega_W 1.0e-7 //relaxation factor for Load
#define omega_P 1.0e-3 //relaxation factor for Pressure
#define Err_Pmin 1.0e-6 //convergence criteria

#define dcrank 2.0 //step size for crank angle

int cen_x;
double c=14.9e-6; //crown height in m
double rw=1.475e-3; //ring width in m
double cr=0.04; //crank radius in m
double conrod=0.1419; //connecting-rod length in m
double bore=0.0889; //bore size in m
double N=2000.0; //engine speed in rpm
double Tau0=2.0e6; //Eyring stress in Pa
double m=0.08; //Pressure coefficient of boundary shear stress
double slope=0.08; // Slope of oil limiting shear stress-pressure relation
double sigma=0.37e-6; //RMS asperity height in m
double  CombP[32][2]={
	{-0.12843e5,-360.0},
	{-0.09476e5,-330.0},
	{-0.09156e5,-300.0},
	{-0.08835e5,-270.0},
	{-0.04827e5,-240.0},
	{-0.03545e5,-210.0},
	{0.01768e5,-180.0},
	{0.07654e5,-150.0},
	{0.15105e5,-120.0},
	{0.32472e5,-90.0},
	{0.85688e5,-60.0},
	{1.52637e5,-45.0},
	{2.97704e5,-30.0},
	{6.87269e5,-15.0},
	{11.4173e5,-5.0},
	{13.62334e5,0.0},
	{15.33501e5,5.0},
	{16.22085e5,12.0},
	{16.05817e5,15.0},
	{12.00915e5,30.0},
	{7.44638e5,45.0},
	{4.59107e5,60.0},
	{2.19067e5,90.0},
	{1.31842e5,120.0},
	{0.83901e5,150.0},
	{0.4535e5,180.0},
	{0.24349e5,210.0},
	{0.11364e5,240.0},
	{0.04631e5,270.0},
	{-0.00529e5,300.0},
	{-0.07537e5,330.0},
	{-0.13484e5,360.0}
	}; //specific combustion pressure

FILE 	*out,*out1,*in;
char name_file[1500], name_dir[1500], name_dir2[1500];

double h0[MAX_x+1], H[MAX_x+1],H_old[MAX_x+1], h_ref;
double P[MAX_x+1], P_old[MAX_x+1], dP[MAX_x+1],term[MAX_x+1], Buf[MAX_x+1];
double u, load;

int main()
{
	cen_x = (MAX_x)/2.0;
	
	int i, loop, loop_w, reversal, first;
	double a, d, circ;
	double x[MAX_x+1];
	double W_discreet, W_set, R_discreet, E1, poi1, E2, poi2, E, b, dx, DX, PH, eta, dh, dt;
	double J0, J1, J2, J3, psi;
    double A0, C0, S_star;
    double Err_P, Err_W, sum_P, diff_P, pmax, hmin;
    double crankangle, Upiston, Upiston0, teta, omega, Dt;
    double z,y,z1,z2,y1,y2;
	double F2, F5_2, lambda, dA;
	double Aa[MAX_x+1], Pa[MAX_x+1], Fb, Fv, TauV[MAX_x+1], TauB[MAX_x+1];
	
	a=rw*0.5;
	circ=2.0*pi*(0.5*bore); //circumference of ring
	eta=6.89e-3; //Viscosity at pressure, P=0 and constant temperature
	S_star=0.0;
    
    W_discreet=0.341e6*circ*rw; //Ring tension
    R_discreet=a;
    
    E1=210.0e9;	poi1=0.30; //Material 1 yield, poisson ratio for material 1 (cylinder wall)
    E2=210.0e9;	poi2=0.30; //Material 2 yield, poisson ratio for material 2 (piston ring)
    
    E=2.0/((1.0-poi1*poi1)/E1+(1.0-poi2*poi2)/E2); //Effective E
    
    d=(4.0*W_discreet)/(pi*E*circ);
    
    b=sqrt(4.0*W_discreet*R_discreet/(pi*E*circ));
    PH = 2.0*W_discreet/(pi*b*circ);
    
   // b=pow(R_discreet*d,0.5); //Contact radius
    //PH=pow((E*W_discreet)/(pi*circ*R_discreet),0.5); //Hertzian pressure
    
    dx=rw/(MAX_x); 
	DX=dx/b; //dimensionless dx
    
    printf("d:%e \n", d);
    printf("Contact radius: %e \n", b);
	printf("Hertzian pressure:%e\n", PH);
	printf("Ring tension:%e\n\n", W_discreet);
	
	//Loop for velocity
	omega=(2.0*pi*N)/60;
	printf("Omega: %e \n\n", omega);
    
    h_ref=1.0e-6;
    for (i=0;i<=MAX_x;i++)
    {
        x[i]=(-rw/2.0)+(i*dx);
        h0[i]=c*pow(x[i]/a,2.0);
        P[i]=0.0;
    }
    
    sprintf(name_file,"%s","./Results/output.txt");
    out = fopen(name_file,"w+");
    fclose(out);
    
    loop_w=0;
    reversal=0;
    first=0;

    Dt=(dcrank*pi/180.0)/omega;
	for (crankangle=-360.0;crankangle<=360.0;crankangle=crankangle+dcrank)
	{
        loop_w++;
        teta=crankangle*pi/180.0;
        Upiston0=Upiston;
		Upiston=cr*omega*(sin(teta)+(0.5*cr*sin(2.0*teta))/conrod);
        
        u=fabs(0.5*Upiston);
        
        if(first==0)
        {
            if(fabs(Upiston)<2.0)
                continue;
            else
                first=1;
        }
        
        else
        {
            if(fabs(Upiston)<0.01) // Reduce to approx 0
                continue;
        }
        
        if(Upiston/Upiston0<0.0)
            reversal=1;
        else
            reversal=0;
    
		y=crankangle;
		
		//Limit for crankangle
		for (i=0;i<32;i++)
		{
			if (y>CombP[i][1])
			{
				y1=CombP[i][1];
				y2=CombP[i+1][1];
	
	            z1=CombP[i][0];
	            z2=CombP[i+1][0];
			}
		}
			
		//Interpolation
		z=((z2-z1)/(y2-y1))*(y-y1)+z1;
		
        if(z<0.0)
            z=0.0;
		//Loop to find set load
		W_set=(0.341e6+z)*circ*rw; //set load = ring tension + combustion pressure
		printf("Set load = %e\n\n",W_set);
		
        for (i=0;i<=MAX_x;i++)
            H_old[i]=H[i];
        
        Err_W=0.0;
        
        //Reversal loop 
        if(reversal==1)
        {
        	printf("\n\t\tEntering reversal loop\n\n");
        	
        	for(i=0;i<=MAX_x;i++)
                {
                    Buf[i]=0.0;
                    Buf[cen_x]=P[cen_x];
                }
            
            for(i=0;i<=MAX_x-1-cen_x;i++)
                {
                    Buf[cen_x-i]=P[cen_x+i];
                    Buf[cen_x+i]=P[cen_x-i];
                }
            
            for(i=0;i<=MAX_x;i++)
                {
                    P[i]=Buf[i];
                    P_old[i]=Buf[i];
                }
            
            sprintf(name_file,"%s%d%s","./Results/",loop_w,"_buffer.txt");
            out = fopen(name_file,"w+");
            
            for (i=1;i<=MAX_x;i++)
            {
                fprintf(out,"%5d %15e %15e",i, x[i],Buf[i]*PH);
                fprintf(out,"\n");
            }
            fclose(out);
            
            printf("Exiting reversal loop\n\n");
        }
            
        do
		{
			loop=1;
        
            psi=12.0*u*eta*pow(R_discreet,2.0)/(PH*pow(b,3.0));
        
            if(hmin>3.0e-6)
                h_ref=h_ref+10.0*Err_W*omega_W;
            else if(hmin>0.5e-6 && hmin<=3.0e-6)
                h_ref=h_ref+Err_W*omega_W;
            else if(hmin>0.1e-6 && hmin<=0.5e-6)
                h_ref=h_ref+0.1*Err_W*omega_W;
            else
                h_ref=h_ref+0.01*Err_W*omega_W;
   
        	//Oil film thickness and squeeze film
            for (i=0;i<=MAX_x;i++)
            {
                H[i]=(h0[i]+h_ref)*R_discreet/(b*b);
                S_star=(H[i]-H_old[i])/(R_discreet/(b*b));
                
                if(S_star>0.0)
                	S_star=0.0;
                    
                term[i]=S_star/(Dt*u);
                     
                    //term[i]=0.0;     
            }
	  	 	
            Err_P=1.0;
            
            do
			{
				for(i=boundary;i<=MAX_x-boundary;i++)
	  		      	P_old[i]=P[i];
		 	           
			   	for(i=boundary;i<=MAX_x-boundary;i++)
	 	 	  	{
			       	A0=(1.0/(2.0*DX*DX))*((pow(H[i+1],3.0)+pow(H[i],3.0))*P[i+1]-(pow(H[i+1],3.0)+2.0*pow(H[i],3.0)+pow(H[i-1],3.0))*P[i]+(pow(H[i],3.0)+pow(H[i-1],3.0))*P[i-1]);     
			       	C0=(1.0/DX)*((1.0-beta)*(H[i+1]*1.0-H[i]*1.0)+beta*(H[i]*1.0-H[i-1]*1.0));
			        	
			       	J0=A0-psi*(C0+R_discreet*term[i]/b);
                    
			       	J1=(1.0/(2.0*DX*DX))*(pow(H[i+1],3.0)+pow(H[i],3.0));
			       	J2=(1.0/(2.0*DX*DX))*(pow(H[i],3.0)+pow(H[i-1],3.0));
			       	J3=-(1.0/(2.0*DX*DX))*(pow(H[i+1],3.0)+2.0*pow(H[i],3.0)+pow(H[i-1],3.0));
			                
			       	dP[i]=(-J0-J1*dP[i+1]-J2*dP[i-1])/J3;
			
                    P[i]=P_old[i]+omega_P*dP[i];
			                
			        if(P[i]<0.0)
			           	P[i]=0.0;
				}
				
				//Convergence criteria
				sum_P=0.0;
		        diff_P=0.0;
		        pmax=0.0;
		        hmin=1.0;
		        
		           for(i=boundary;i<=MAX_x-boundary;i++)
		            {
		                if(P[i]>0.0)
		                {
			                sum_P+=fabs(P[i]);
				            diff_P+=fabs(P[i]-P_old[i]);
			            }
			            if (P[i]>pmax)
			                pmax = P[i];
	
						if(H[i]*b*b/R_discreet<hmin)
		                 	hmin = H[i]*b*b/R_discreet;
			        }
			        
			    Err_P=diff_P/sum_P;
			           
			        sprintf(name_file,"%s%d%s","./Results/",loop_w,"result.txt");
			        out = fopen(name_file,"w+");
			        for (i=0;i<=MAX_x;i++)
			        {
			          	fprintf(out,"%5d %15e %15e %15e %15e %15e",i,x[i],H[i]*b*b/R_discreet,P[i]*PH, term[i],H_old[i]*b*b/R_discreet);
			            fprintf(out,"\n");
			        }
			        fclose(out);
			    
			        if (loop%1000==0)
			       		printf("\t%d%15e\n",loop,Err_P); //must become smaller, must reach Pmin
						         
			        loop++;
			}
            while (Err_P>Err_Pmin || loop<=100);
				    
			// This is to calculate the load based on the pressure
			load=0.0;
			for(i=boundary;i<=MAX_x-boundary;i++)
				load+=(P[i]*PH)*dx*circ;
				
				Err_W=(load-W_set)/W_set;
				printf("\n\t\tVelocity=%e\t h_ref=%e\tErr_W=%e\tload=%e\tcrankangle=%15e\n\n",Upiston,h_ref,(Err_W),load,crankangle);
				
		}
		while(fabs(Err_W)>1.0e-3);
		
		//Rough surface contact model
		
        Fb=0.0;
        Fv=0.0;
        for(i=boundary;i<=MAX_x-boundary;i++)
        {
            lambda=(H[i]*b*b/R_discreet)/sigma;
            
			F5_2=-0.0046*pow(lambda,5.0)+0.0574*pow(lambda,4.0)-0.2958*pow(lambda,3.0)+0.7844*pow(lambda,2.0)-1.0776*lambda+0.6167; //F5/2
            F2=-0.0018*pow(lambda,5.0)+0.0281*pow(lambda,4.0)-0.1728*pow(lambda,3.0)+0.5258*pow(lambda,2.0)-0.8043*lambda+0.5003; //F2	
	
            if (F5_2<=0.0)
                F5_2=0.0;
		
            if (F2<=0.0)
                F2=0.0;
		
            Aa[i]=0.0298*F2*dx*circ; //asperity area (m2)
            Pa[i]=0.000227*F5_2*E*dx*circ; //asperity load (N)
            
            if(Aa[i]>0)
                TauB[i]=(Tau0*Aa[i]+m*Pa[i])/(dx*circ); //shear for boundary friction
            else
                TauB[i]=0.0;
		
            TauV[i]=eta*fabs(Upiston)/(H[i]*b*b/R_discreet); //shear for viscous friction
            
            if (TauV[i]>=Tau0)
                TauV[i]=Tau0+slope*P[i];
        
            dA=(dx*circ)-Aa[i];
		
			Fb+=TauB[i]*(dx*circ); //Boundary friction force
			Fv+=TauV[i]*dA; //Viscous friction force
		}
        
        
        Fb=Fb*Upiston/fabs(Upiston);
        Fv=Fv*Upiston/fabs(Upiston);
        
        sprintf(name_file,"%s%d%s","./Shear/",loop_w,"result.txt");
        out = fopen(name_file,"w+");
        for (i=0;i<=MAX_x;i++)
        {
            fprintf(out,"%5d %15e %15e %15e",i,x[i],TauV[i],TauB[i]);
            fprintf(out,"\n");
        }
        fclose(out);
		
		sprintf(name_file,"%s","./Results/output.txt");
		out=fopen(name_file,"a");
		fprintf(out,"%d %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e \n",loop_w,crankangle,Upiston,W_set,z,load,hmin,pmax*PH,term[cen_x],Fb,Fv, Fb+Fv);
		fclose(out);
	}

}
