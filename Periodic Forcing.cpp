
#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>
#include <stdlib.h>    
#include <stdio.h>


//R_A_N_2_P_A_R_A_M_E_T_E_R_S___________________________________________________

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 3.14159265359

using namespace std;

#define T 1024*4

#define J 1000///discretization for theoretical frequency calculations
#define Q 1//Discretization for parameter scanning
#define P 1//Trials
#define dt 0.1
#define Dt 0.001
#define N 100
#define alpha 1.0
#define gain 300
#define h 0.0
#define sparseness 0.
#define g -10
#define delay 25
#define R_c -1.375
#define D_max 0.1

// PROBLEM alpha^-1 *dx/dt = -x+R_1*x(t-/tau) 


double RealEigenproblem(double rel, double iml,double R_1, double tau);
double Heaviside(double u);
double F_corrected(double u, double A, double w);
double F_Axel(double u, double A, double w);
double ImaginaryEigenproblem(double rel, double iml,double R_1, double tau);
float ran2(long *idum);//random number generator - initial cond.
void shuffle(double Array[], int size);//shuffle an array with double precision entries
void four1(double data[], unsigned long nn, int isign);// FFT NR routine
double signum(double u);
double f(double u);
double f_prime(double u);
double f_2prime(double u);
double f_3prime(double u);
double f_4prime(double u);
double f_5prime(double u);
double f_7prime(double u);
double f_9prime(double u);
double f_11prime(double u);
double u[N][T];
double u_spike[N][T];
double 	u_spatial_ave[T]; 
double u_ave[T];
double u_lin[T];
double u_approx[T];
double u_ave_spike[T];
double u_mean[T];
double u_fixed_point[T];
double nu[T];
double Dummy_1[T];
double Dummy_2[T];
double Dummy_3[T];
double Dummy_4[T];
double v[N][T];
double F[N][T];
double F2[N][T];
double F3[N][T];
double X[N][T];
double XI[N][T];
double W[N][N];
double PSD_1[(int)T];//Power Spectral Density 
double PSD_2[(int)T];//Power Spectral Density 
double PSD_3[(int)T];//Power Spectral Density 
double PSD_4[(int)T];//Power Spectral Density 
double Freq[T/2];//Array of output frequencies for PSD display
double D[Q];
double I_o[Q];
double freq[Q];
double PERIODIC[N][T];
double Peak_frequency_1[Q];
double Peak_frequency_2[Q];
double Peak_frequency_3[Q];
double Peak_frequency_4[Q];
double Theoretical_frequency[Q];
double Susceptibility[Q];
double Susceptibility_theo[Q];
double Fixed_point[Q];
double Fixed_point_theo[Q];

int main()
{
	double tau= delay*dt;
	double mean_w=0;
	for (int n=0;n<N;n++)
	{
			
			for (int l=0;l<N;l++)
			{
				long seed_w=(long) n*l+l+n+34;
				long seed=(long)  32+n*l+n+l+12;
				
						if(ran2(&seed)<sparseness)
						{
								W[n][l]=0;
							//	W[l][n]=W[n][l];
						}
						//else{W[n][l] =  -g*ran2(&seed_w);	/*W[l][n]=W[n][l];*/	}
					   else{W[n][l]=g+2*fabs(g)*sqrt(-2*log(ran2(&seed_w)))*cos(2*3.14159*ran2(&seed));}
				
				mean_w=mean_w+1/(((double) N)*((double)N))*W[n][l];
			}
		
	}
//	cout<<mean_w<<endl;



for (int q=0;q<Q;q++)
{
	cout<<"WARNING: includes periodic forcing!"<<endl;
	D[q] =0.0;//+ D_max*q/((double)Q-1);
	I_o[q]=1.0;//+1*q/((double)Q);
	freq[q]=100;
	cout<<freq[q]<<endl;
	double omega = freq[q];//(2*Pi)*freq[q]*Dt/dt;
	cout<<I_o[q]/freq[q]/*I_o[q]*I_o[q]*Pi/(((omega*omega)+1)*omega)*/<<endl;	

	



for (int p=0;p<P;p++)
{

	

				double var_nu=0;
				for (int t=0;t<T;t++)
				{
					for (int n=0;n<N;n++)
					{
					
						 long d=rand()%105;      
			             long noisyseed1=(long)21*t+d*n+q*p+90*p+45*q;
			             long noisyseed2=(long)69*t+11*n+d*t+q*p+901*p+415*q; 
			             XI[n][t] =  sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
			             PERIODIC[n][t] = I_o[q]*sin(2*Pi*freq[q]*t*Dt);
			        }
					 nu[t+1]=nu[t]+dt*((-nu[t])+I_o[q]*PERIODIC[0][t])+sqrt(2*D[q]*dt)*XI[0][t];
					var_nu=var_nu+1/((double)T)*(nu[t+1]*nu[t+1]);	  
						  
				}
			//	cout<<"var nu:"<<var_nu<<endl;
					
				for (int t=tau;t<T;t++)
				{
						
						
						for (int i=0;i<N;i++)
						{
							
							u[i][t+1]= u[i][t]+dt*(-u[i][t]+F[i][t-delay]+PERIODIC[i][t]);
							v[i][t+1]= v[i][t]+dt*(-v[i][t]+F2[i][t-delay]);
							u_spike[i][t+1]=u_spike[i][t]+dt*(-u_spike[i][t]+F3[i][t-delay]+PERIODIC[i][t]);
						
							F[i][t]=0;
							F2[i][t]=0;
							F3[i][t]=0;
							for (int j=0;j<N;j++)
							{
								F[i][t] = F[i][t] + 1/((double)N)*W[i][j]*f(u[j][t]);
							}
							for (int j=0;j<N;j++)
							{
								F2[i][t] = F2[i][t] + 1/((double)N)*W[i][j]*f(v[j][t]);
							}
							for (int j=0;j<N;j++)
							{
								F3[i][t] = F3[i][t] + 1/((double)N)*W[i][j]*X[j][t];
							}
							
							long seed =rand()%101;  
							long seed_1 = (long) 3+seed*i+11*q+45*p*i+i+t;                                                  
							double p =ran2(&seed_1);
							double p_fire =  1-exp(-f(u_spike[i][t+1])*dt);                         
							if (p<p_fire)//spike occurs
							{                   
							                 X[i][t+1] =1/dt;                                  
							}
							else//no spike
							{
									          X[i][t+1]=0; 
							}
							
							u_spatial_ave[t+1] = 	u_spatial_ave[t+1]+1/((double)N)*u[i][t+1];
							
							
						}
						u_ave_spike[t]=0;
						for (int i=0;i<N;i++)
						{
							u_ave_spike[t]=u_ave_spike[t]+1/((double)N)*u_spike[i][t+1];
						}
						
						double A1 = 0;
						double B = 0;
						
						for (int s=0;s<T;s++)
						{
							
							A1 = A1 + 1/((double) T)*f(u_ave[t-delay]+nu[s]);
							B = B + 1/((double) T)*f(u_fixed_point[t]+nu[s]);
						
						}
						//cout<<C<<endl;
						
						
						u_ave[t+1]= u_ave[t]+dt*(-u_ave[t]+mean_w*A1);
					
				
			
						
					
						u_mean[t+1]= u_mean[t]+dt*(-u_mean[t]+mean_w*(F_Axel(u_mean[t-delay], I_o[q],freq[q]/(2*Pi)))+1*PERIODIC[0][t]);
						//u_mean[t+1]= u_mean[t]+dt*(-u_mean[t]+mean_w*(f(u_mean[t-delay]))+PERIODIC[0][t]);
					
						u_fixed_point[t+1] = u_fixed_point[t]+dt*(-u_fixed_point[t]+mean_w*B);
						
				}
					
				
				
				double R;
				
				double c_min=-I_o[q]/freq[q]+0.00001;
				int Nc =1000;
				double c_max=I_o[q]/freq[q]-0.00001;
				double dc = fabs(c_min-c_max)/((double) Nc);
				double C=0;
				double E=0;
				double F=0;
				double G=0;
				
				
					for (int k=0;k<Nc;k++)
					{
										 double c = c_min+dc*k; 
										 C = C + dc*mean_w*f_prime(u_fixed_point[T-10]+c)*1/sqrt((I_o[q]/freq[q]*I_o[q]/freq[q])-c*c)*1/Pi;
									    
									
									
					}
			
			
				R=C;
			
				Fixed_point[q]=Fixed_point[q]+1/((double)P)*u_fixed_point[T-10];
				Fixed_point_theo[q]=0;//-sqrt(Pi)*sqrt(D[q])*mean_w/(sqrt(2)*g-2*sqrt(Pi)*sqrt(D[q]));
				//	Susceptibility[q]=Susceptibility[q]+1/((double)P)*R;
					Susceptibility[q]=Susceptibility[q]+1/((double)P)*(C);
			//	Susceptibility_theo[q]=mean_w*(Heaviside(Fixed_point[q]+l)/(sqrt(0.1e-1-(Fixed_point[q]*Fixed_point[q])+l^2)*Pi)-Heaviside(mu-l)/(sqrt(0.1e-1-mu^2+l^2)*Pi))
			//	cout<<"R:"<<Susceptibility[q]<<endl;
			//	cout<<"fp:"<<Fixed_point[q]<<endl;
				
					for (int t=tau;t<T;t++)
				{
						
					u_lin[t+1]=u_lin[t]+dt*(-u_lin[t]+R*u_lin[t-delay]);
						
				}
				
			    	//Compute mean EEG
											//Compute mean EEG
							double mean_1=0;
							double mean_2=0;
							double mean_3=0;
							double mean_4=0;
							for (int t=0;t<(int)T;t++)
							{
														mean_1  = mean_1+1/((double) T)*u[(int)N/2][t];
														mean_2  = mean_2+1/((double) T)*v[(int)N/2][t];
														mean_3  = mean_3+1/((double) T)*u_mean[t];
														mean_4  = mean_4+1/((double) T)*u_ave_spike[t];
							}
										
														
											//Compute mean corrected EEG signal
											for (int t=0;t<(int)T;t++)
											{
														Dummy_1[t]  = u[(int)N/2][t]-mean_1;
														Dummy_2[t]  = v[(int)N/2][t]-mean_2;
														Dummy_3[t]  = u_mean[t]-mean_3;
														Dummy_4[t]  = u_ave_spike[t]-mean_4;
											}
										
											//fourier transforms and power spectral denstities
											for (int k=0;k<(int)T;k++)
											{
														PSD_1[k] = 0;	
														PSD_2[k] = 0;	
														PSD_3[k] = 0;	
														PSD_4[k] = 0;	 		 
										    }
																			
											for (int k=0;k<(int)T;k++)
											{
														PSD_1[k] = Dummy_1[k];
														PSD_2[k] = Dummy_2[k];
														PSD_3[k] = Dummy_3[k];
														PSD_4[k] = Dummy_4[k];
															 		 
											}
											
											unsigned long nn=T/2;
											four1(PSD_1-1, nn,1);
										
											 
											for (int k=0;k<T/2;k++)
											{
													PSD_1[k] = 1/((double)T*T)*(fabs(PSD_1[k])*fabs(PSD_1[k])+fabs(PSD_1[(int)T-k])*fabs(PSD_1[(int)T-k]));
											}
											
										
											four1(PSD_2-1, nn,1);
										
											 
											for (int k=0;k<T/2;k++)
											{
													PSD_2[k] = 1/((double)T*T)*(fabs(PSD_2[k])*fabs(PSD_2[k])+fabs(PSD_2[(int)T-k])*fabs(PSD_2[(int)T-k]));
											}
											
												four1(PSD_3-1, nn,1);
										
											 
											for (int k=0;k<T/2;k++)
											{
													PSD_3[k] = 1/((double)T*T)*(fabs(PSD_3[k])*fabs(PSD_3[k])+fabs(PSD_3[(int)T-k])*fabs(PSD_3[(int)T-k]));
											}
											
											four1(PSD_4-1, nn,1);
											
											for (int k=0;k<T/2;k++)
											{
													PSD_4[k] = 1/((double)T*T)*(fabs(PSD_4[k])*fabs(PSD_4[k])+fabs(PSD_4[(int)T-k])*fabs(PSD_4[(int)T-k]));
											}
											
											double max_f1=0;
											double max_power1=0;
											double max_f2=0;
											double max_power2=0;
											double max_f3=0;
											double max_power3=0;
											double max_f4=0;
											double max_power4=0;
											for(int k =2*2*(T)*Dt; k<2*50*(T)*Dt;k++)
											{
													if (	PSD_1[k] >max_power1)
													{
														max_power1 = PSD_1[k];
														max_f1=k*1/((double) 2* (T)*Dt);
													}	
											}
											for(int k =2*2*(T)*Dt; k<2*50*(T)*Dt;k++)
											{
													if (	PSD_2[k] >max_power2)
													{
														max_power2 = PSD_2[k];
														max_f2=k*1/((double) 2* (T)*Dt);
													}	
											}
											for(int k =2*2*(T)*Dt; k<2*50*(T)*Dt;k++)
											{
													if (	PSD_3[k] >max_power3)
													{
														max_power3 = PSD_3[k];
														max_f3=k*1/((double) 2* (T)*Dt);
													}	
											}
											for(int k =2*2*(T)*Dt; k<2*50*(T)*Dt;k++)
											{
													if (	PSD_4[k] >max_power4)
													{
														max_power4 = PSD_4[k];
														max_f4=k*1/((double) 2* (T)*Dt);
													}	
											}
											Peak_frequency_1[q] = Peak_frequency_1[q]+1/((double) P)*max_f1;
											Peak_frequency_2[q] = Peak_frequency_2[q]+1/((double) P)*max_f2;
											Peak_frequency_3[q] = Peak_frequency_3[q]+1/((double) P)*max_f3;
											Peak_frequency_4[q] = Peak_frequency_4[q]+1/((double) P)*max_f4;
											Theoretical_frequency[q] = -(dt/Dt)/(2*Pi)*(-Pi+acos((-1)/R_c))/(delay*dt);;//-2*Pi*R*sqrt(1-(-1/R)*(-1/R));//
											cout<<Peak_frequency_1[q]<<"	"<<Peak_frequency_2[q]<<endl;
	}//p loop
							
	}//q loop
							
	double H=0;	
	int No_of_points=100;
	double U_max=0.5;
	double U_min=-0.5;
	double U[No_of_points];
	double F[No_of_points];
	double F_approx[No_of_points];
	
	int Nc =100000;
	
	double Input_Level=0.1;
	double fre =2*Pi*50/100;
	double c_min=-Input_Level/fre+0.001*Input_Level;

	double c_max=Input_Level/fre-0.001*Input_Level;
	double dc = fabs(c_min-c_max)/((double) Nc);
	for (int k1=0;k1<100;k1++)
	{
		U[k1] = U_min+fabs(U_max-U_min)/((double)No_of_points)*k1;
		H=0;
		for (int k2=0;k2<Nc;k2++)
		{
							 double c = c_min+dc*k2; 
						//	 cout<<c<<endl;
							 H = H + dc*f(U[k1]+c)*1/(Pi*sqrt(Input_Level/fre*Input_Level/fre-c*c));
							 
						
						
		}
		cout<<H<<endl;
		F[k1] = H;
		//F_approx[k1]=0.5+(0.5)*erf((0.5)*sqrt(2)*U[k1]/sqrt(Noise_Level));
	}		
			
			
				
					for(int k=0;k<(int) T/2;k++)
					
					{
								Freq[k] = k*1/((double) 2* (T)*Dt);
								//cout<<k*1/((double) 2* (T/2)*Dt)<<endl;
					}


    ofstream outfile;
    
     outfile.open("SS  -  PSD.txt", ios::out);
    for(int k=0;k<(int) 30*(2*T)*Dt;k++)
    {
         	
         		  			outfile<<Freq[k]<<"	"<<PSD_1[k]<<"	"<<PSD_2[k]<<endl;  	
    }  
    outfile.close(); 
    
      outfile.open("SS  -  Connectivity.txt", ios::out);
    for(int i=0;i<N;i++)
    {
    	for (int j=0;j<N;j++)
         {
			
         		  			outfile<<i<<"	"<<j<<"	"<<W[i][j]<<endl;  	
        }
    }  
    outfile.close(); 
    
     outfile.open("SS  -  Noise induced Linearization.txt", ios::out);
    for(int k=0;k<No_of_points;k++)
    {
         	
         		  			outfile<<U[k]<<"	"<<F[k]<<"	"<<F_approx[k]<<endl;  	
    }  
    outfile.close(); 
    
     outfile.open("SS - Sliced view.txt", ios::out);
     for(int i=0;i<T;i++)
     {
                           outfile<<i*Dt<<"	"<<u[(int)N/2][i]<<"	"<<v[(int)N/2][i]<<"	"<<	u_spatial_ave[i]<<"	"<<u_mean[i]<<"	"<<u_ave_spike[i]<<"	"<<I_o[0]*PERIODIC[(int)N/2][i]<<endl;
						  
						  
	}
           
     outfile.close(); 
     
      
     outfile.open("SS - misceleanous view.txt", ios::out);
     for(int i=0;i<T;i++)
     {
                           outfile<<i*Dt<<"	"<<u[(int)N/2][i]<<"	"<<u[(int)N/2+15][i]<<"	"<<u[(int)N/2+9][i]<<"	"<<u[(int)N/2-15][i]<<"	"<<u[(int)N/2+10][i]<<endl;
						  
						  
	}
           
     outfile.close(); 
     
     outfile.open("SS - Spikes.txt", ios::out);
     for(int i=0;i<T;i++)
     {
     	for (int n=0;n<N;n++)
     	{
		 
     			if (X[n][i])
		     	{
		     		 outfile<<i*Dt<<"	"<<n<<endl;
				}
                          
		}
						  
	}
           
     outfile.close(); 
        
     outfile.open("SS - Peak frequency.txt", ios::out);
     for(int q=0;q<Q;q++)
     {
                           outfile<<I_o[q]/freq[q]<<"	"<<Peak_frequency_1[q]<<"	"<<Peak_frequency_2[q]<<"	"<<Peak_frequency_3[q]<<"	"<<Peak_frequency_4[q]<<"	"<<Theoretical_frequency[q]<<endl;
						  
						  
	}
	
 	outfile.close(); 
 	
 	   outfile.open("SS - Corrected F with Periodic Forcing.txt", ios::out);
     for(int q=0;q<1000;q++)
     {
     	double input=-0.5+1*q/((double)1000);
                           outfile<<input<<"	"<<F_corrected(input, 0.01,2*Pi*50/(100))<<endl;
						  
	}
 	outfile.close(); 
 	
 	 outfile.open("SS - R.txt", ios::out);
     for(int q=0;q<Q;q++)
     {
                           outfile<<I_o[q]<<"	"<<Susceptibility[q]<<"	"<<Susceptibility_theo[q]<<endl;
						  
						  
	}
 	outfile.close(); 
 	
 	 outfile.open("SS - Fixed Point.txt", ios::out);
     for(int q=0;q<Q;q++)
     {
                           outfile<<D[q]<<"	"<<Fixed_point[q]<<"	"<<Fixed_point_theo[q]<<endl;
						  
						  
	}
 	outfile.close(); 
    system("pause");
    cout<<endl;
      
return 0;    
}

double RealEigenproblem(double rel, double iml, double R_1, double tau)
{
       double output=rel*(1/alpha)+1+R_1*exp(-rel*tau)*cos(iml*tau);
       
       return output;
       }
       
double ImaginaryEigenproblem(double rel, double iml,double R_1, double tau)
{
       double output=iml*(1/alpha)-R_1*exp(-rel*tau)*sin(iml*tau);
       
       return output;
       }


double signum(double u)
{
	double output;
	if (u>0)
	{
		output=1;
		
	}
	else{output=-1;}
	return output;
}

double f(double u)
{
	double output;
	output = 1/(1+exp(-gain*(u-h)));
	//output = tanh(u);//1/(1+exp(-gain*(u-h)));
	return output;
	
}

double f_prime(double u)
{
	double output;
	output = pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * gain * exp(-gain * (u - h));
	//output = 1-tanh(u)*tanh(u);//pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * gain * exp(-gain * (u - h));
	return output;
	
}

double f_2prime(double u)
{
	double output;
	output = 0.2e1 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * gain * gain * pow(exp(-gain * (u - h)), 0.2e1) - pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * gain * gain * exp(-gain * (u - h));
	//output = -0.2e1 * tanh(u) * (0.1e1 - pow(tanh(u), 0.2e1)); //0.2e1 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * gain * gain * pow(exp(-gain * (u - h)), 0.2e1) - pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * gain * gain * exp(-gain * (u - h));
	return output;
	
}
double f_3prime(double u)
{
	double output;
	output =  0.6e1 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.3e1) * pow(exp(-gain * (u - h)), 0.3e1) - 0.6e1 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.3e1) * pow(exp(-gain * (u - h)), 0.2e1) + pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.3e1) * exp(-gain * (u - h));
	//output = -0.2e1 * pow(0.1e1 - pow(tanh(u), 0.2e1), 0.2e1) + 0.4e1 * pow(tanh(u), 0.2e1) * (0.1e1 - pow(tanh(u), 0.2e1));// 0.6e1 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.3e1) * pow(exp(-gain * (u - h)), 0.3e1) - 0.6e1 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.3e1) * pow(exp(-gain * (u - h)), 0.2e1) + pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.3e1) * exp(-gain * (u - h));
	return output;
	
}
double f_5prime(double u)
{
	double output;
	output = 0.120e3 * pow(0.1e1 + exp(-gain * (u - h)), -0.6e1) * pow(gain, 0.5e1) * pow(exp(-gain * (u - h)), 0.5e1) - 0.240e3 * pow(0.1e1 + exp(-gain * (u - h)), -0.5e1) * pow(gain, 0.5e1) * pow(exp(-gain * (u - h)), 0.4e1) + 0.150e3 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.5e1) * pow(exp(-gain * (u - h)), 0.3e1) - 0.30e2 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.5e1) * pow(exp(-gain * (u - h)), 0.2e1) + pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.5e1) * exp(-gain * (u - h));
	return output;
	
}
double f_4prime(double u)
{
	double output;
	output = 0.24e2 * pow(0.1e1 + exp(-gain * (u - h)), -0.5e1) * pow(gain, 0.4e1) * pow(exp(-gain * (u - h)), 0.4e1) - 0.36e2 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.4e1) * pow(exp(-gain * (u - h)), 0.3e1) + 0.14e2 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.4e1) * pow(exp(-gain * (u - h)), 0.2e1) - pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.4e1) * exp(-gain * (u - h));

		return output;
	
}



double f_7prime(double u)
{
	double output;
	output = 0.5040e4 * pow(0.1e1 + exp(-gain * (u - h)), -0.8e1) * pow(gain, 0.7e1) * pow(exp(-gain * (u - h)), 0.7e1) - 0.15120e5 * pow(0.1e1 + exp(-gain * (u - h)), -0.7e1) * pow(gain, 0.7e1) * pow(exp(-gain * (u - h)), 0.6e1) + 0.16800e5 * pow(0.1e1 + exp(-gain * (u - h)), -0.6e1) * pow(gain, 0.7e1) * pow(exp(-gain * (u - h)), 0.5e1) - 0.8400e4 * pow(0.1e1 + exp(-gain * (u - h)), -0.5e1) * pow(gain, 0.7e1) * pow(exp(-gain * (u - h)), 0.4e1) + 0.1806e4 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.7e1) * pow(exp(-gain * (u - h)), 0.3e1) - 0.126e3 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.7e1) * pow(exp(-gain * (u - h)), 0.2e1) + pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.7e1) * exp(-gain * (u - h));
	return output;
	
}
double f_9prime(double u)
{
	double output;
	output = 0.362880e6 * pow(0.1e1 + exp(-gain * (u - h)), -0.10e2) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.9e1) - 0.1451520e7 * pow(0.1e1 + exp(-gain * (u - h)), -0.9e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.8e1) + 0.2328480e7 * pow(0.1e1 + exp(-gain * (u - h)), -0.8e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.7e1) - 0.1905120e7 * pow(0.1e1 + exp(-gain * (u - h)), -0.7e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.6e1) + 0.834120e6 * pow(0.1e1 + exp(-gain * (u - h)), -0.6e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.5e1) - 0.186480e6 * pow(0.1e1 + exp(-gain * (u - h)), -0.5e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.4e1) + 0.18150e5 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.3e1) - 0.510e3 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.9e1) * pow(exp(-gain * (u - h)), 0.2e1) + pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.9e1) * exp(-gain * (u - h));
	return output;
	
}
double f_11prime(double u)
{
	double output;
	output =  0.39916800e8 * pow(0.1e1 + exp(-gain * (u - h)), -0.12e2) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.11e2) - 0.199584000e9 * pow(0.1e1 + exp(-gain * (u - h)), -0.11e2) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.10e2) + 0.419126400e9 * pow(0.1e1 + exp(-gain * (u - h)), -0.10e2) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.9e1) - 0.479001600e9 * pow(0.1e1 + exp(-gain * (u - h)), -0.9e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.8e1) + 0.322494480e9 * pow(0.1e1 + exp(-gain * (u - h)), -0.8e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.7e1) - 0.129230640e9 * pow(0.1e1 + exp(-gain * (u - h)), -0.7e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.6e1) + 0.29607600e8 * pow(0.1e1 + exp(-gain * (u - h)), -0.6e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.5e1) - 0.3498000e7 * pow(0.1e1 + exp(-gain * (u - h)), -0.5e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.4e1) + 0.171006e6 * pow(0.1e1 + exp(-gain * (u - h)), -0.4e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.3e1) - 0.2046e4 * pow(0.1e1 + exp(-gain * (u - h)), -0.3e1) * pow(gain, 0.11e2) * pow(exp(-gain * (u - h)), 0.2e1) + pow(0.1e1 + exp(-gain * (u - h)), -0.2e1) * pow(gain, 0.11e2) * exp(-gain * (u - h));
	return output;
	
}
float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}


/******************************************************************************/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

double Heaviside(double u)
{
	
	double output;
	if (u<0)
	{
		output=0;
	}
	else
	{
		output=1;
		}	
		return output;
   }   

double F_corrected(double u, double A, double w)
{
	
	double output;
	


	if (u>-A/w&&u<A/w)
	{
		output=(0.5)+asin(w*u/A)/Pi;
	}
	else
	{
		if (u<0)
		{
			output=0;
		}
		else
		{
			output=1;
		}
	}
	
	return output;
	
	
}

double F_Axel(double u, double A, double w)
{
	

              double output;
              double mu=A/w;
              if (u>-A/w&&u<A/w)
              {
                            output=0.5*((2/Pi)*asin(u/mu)+1);
              }
              else
              {
                            if (u<-mu)
                            {
                                           output=0;
                            }
                            else
                            {
                                           output=1;
                            }
              }
              return output;
}

	
	

