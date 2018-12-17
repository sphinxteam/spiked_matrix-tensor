/*****************************************************
The code solves numerically the  Crisanti-Horner-Sommers-Cugliandolo-Kurchan equations 
with a fixed grid, refer to [Sarao Mannelli, S., Biroli, G., Cammarota, C., Krzakala, F., Urbani, P., & Zdeborova, L. (2018)] for details.
All the integrals are evalutated using the simple 2 points Newton-Cotes formula (coefficient 1/2,1/2).
All the derivatives are solved using the Euler method (f(t+dt)=f(t)+dt*(the rest)).

written: 10/12/18 by Stefano Sarao Mannelli and Pierfrancesco Urbani.
******************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

// This structure contains the physical parameters of the system.
typedef struct _psys {
	float p,Delta2,Deltap,Tg,alpha,Cbar0;	// system parameters
	float W22_0,Wp2_0,r2_0,rp_0;			// auxiliary variables
	float W22,Wp2,r2,rp;					// auxiliary variables
	float *T2, *Tp;							// mismatching parameters
} psys;


//	The function computes the equation for mu.
float compute_mu(float*mu, float*Cbar, float**C,float** R, int t, psys *w,float h);
//	The function computes the equation for Cbar. 
float compute_Cbar(float*mu, float*Cbar, float**C,float** R, int t, psys *w,float h);
//	The function computes the equation for C.
float compute_C(int t, int l,float*mu, float*Cbar, float**C,float** R, psys *w,float h);
//	The function computes the equation for R.
float compute_R(int t, int l,float*mu, float*Cbar, float**C,float** R, psys *w,float h);

//	The function produces the output file concerning correlation and response function.
float print(float** C,float** R,float* mu,float *Cbar,int Nmax,float h,FILE* correlation_file,FILE* response_file,FILE* mu_file,FILE* time_file);
//	The function produces the output file concerning the Lagnange multiplier.
float print_mu(float** C,float** R,float* mu,float *Cbar,int Nmax,float h,FILE* mu_file);

//	The function computes the effective noise of the tensor channel.
void dumping_policy(psys *w,float C, float tau, int Nmax, float h, bool protocol);

//	The function propagates the initial values of C,R,Cbar in time up to a given time.
int main (int argc, char *argv[]) {	
	psys w;
	float h, tmax, Cbar0;
	float **C, **R, *Cbar, *mu;
	int Nmax,i,j,t,l;
	char suffix[100], directory[100]="";
	char C_file[100], R_file[100], muname_file[100], T_file[100];
	FILE* correlation_file;	FILE* response_file; FILE* mu_file; FILE*time_file;
	
	// Assign the parameters of the system.
	w.p=3.;
	w.Delta2=.1;			// the variance of the matrix channel
	w.Deltap=.2;			// the variance of the tensor channel
	w.Cbar0 = 0.0001;		// initial value of Cbar
	Nmax=1000;				// number of points in the discretization of the grid
	tmax=10.;				// maximum simulation time

	h = tmax/((float)Nmax)/w.Delta2; // time step size 
	// auxialiary varaibles
	w.W22_0=w.Delta2; w.Wp2_0=w.Delta2*w.Delta2/w.Deltap*2./w.p;
	w.r2_0=1.; w.rp_0=w.Delta2/w.Deltap*2./w.p;
	w.W22 = w.W22_0; w.Wp2 = w.Wp2_0;
	w.r2 = w.r2_0; w.rp = w.rp_0;

	// Use protocol? true - false.
	bool protocol=true;
	float C_policy = 100., tau_policy = 70.;


	fprintf(stderr, "Dynamics for (2+%.0f)-spin planted with Delta2 %.4f and Deltap %.4f (alpha %.4f)\n",w.p,w.Delta2,w.Deltap,w.alpha);
	fprintf(stderr, "System initialized to %.2e magnetization\n",w.Cbar0);
	fprintf(stderr, "Maximum simulation time %.0f using %d steps of size %.2e\n",tmax,Nmax,h);

	sprintf(suffix, "P%.0f_Policy%.0f_%.0f_Delta2-%.4f_DeltaP-%.4f_Cbar-%.2e_t-%.2f_Nmax-%d",w.p,C_policy,tau_policy,w.Delta2,w.Deltap,w.Cbar0,tmax,Nmax);
	sprintf(C_file,"%sC_%s.txt",directory,suffix);
	sprintf(R_file,"%sR_%s.txt",directory,suffix);
	sprintf(muname_file,"%smu_%s.txt",directory,suffix);
	sprintf(T_file,"%sT_%s.txt",directory,suffix);

	// memory allocation.
	correlation_file=fopen(C_file,"w");
	response_file=fopen(R_file,"w");
	time_file=fopen(T_file,"w");
	mu_file=fopen(muname_file,"w");
	C=(float**)malloc((Nmax+1)*sizeof(float*));
	R=(float**)malloc((Nmax+1)*sizeof(float*));
	Cbar=(float*)malloc((Nmax+1)*sizeof(float));
	mu=(float*)malloc((Nmax+1)*sizeof(float));
	w.T2=(float*)malloc((Nmax+1)*sizeof(float));
	w.Tp=(float*)malloc((Nmax+1)*sizeof(float));
	for(i=0; i<=Nmax;i++){  
		C[i]=(float*)malloc((i+1)*sizeof(float));
		R[i]=(float*)malloc((i+1)*sizeof(float));
	}
	
	// inizialization.
	C[0][0]=1;
	R[0][0]=0;
	Cbar[0]=w.Cbar0;
	
	// define protocol.
	dumping_policy(&w,C_policy,tau_policy,Nmax,h,protocol);

	// integration loop.
	for(i=0; i<=Nmax-1; i++){	
		t=i;
		compute_mu(mu, Cbar, C, R, t, &w, h);

		C[t+1][t+1]=1;
		R[t+1][t+1]=0;

		Cbar[t+1]=compute_Cbar(mu, Cbar, C, R, t, &w, h);
		// compute the two times observables.
		for(l=0;l<=t && t<Nmax;l++){
			C[t+1][l]=compute_C(t,l,mu,Cbar,C,R,&w,h);
			R[t+1][l]=compute_R(t,l,mu,Cbar,C,R,&w,h);
		}

		// print progress.
		if(i%100==0){
			fprintf(stderr, "Iteration %d of %d\n",i,Nmax);
		}
		R[t+1][t]=1;
	}
	

	print(C,R,mu,Cbar,Nmax,h*w.Delta2,correlation_file,response_file,mu_file,time_file);
	print_mu(C,R,mu,Cbar,Nmax,h*w.Delta2,mu_file);

	fclose(correlation_file); fclose(response_file);
	fclose(time_file); fclose(mu_file);
	free(C); free(R); free(mu);
		
	return 0;	
}

float compute_mu(float*mu, float*Cbar, float**C,float** R, int t, psys *w,float h){
	int l;
	float bg;
	float auxp,aux2;
	
	aux2  = C[t][0]*R[t][0]/w->T2[0];
	auxp  = 0.5*w->p*w->p*0.5*pow(C[t][0],w->p-1)*R[t][0]/w->Tp[0];
	for(l=1;l<=t-1;l++){
		aux2  += 2.*C[t][l]*R[t][l]/w->T2[l];
		auxp  += w->p*w->p*0.5*pow(C[t][l],w->p-1)*R[t][l]/w->Tp[l];
	}
	aux2 += C[t][t]*R[t][t]/w->T2[t];
	auxp += 0.5*w->p*w->p*0.5*pow(C[t][t],w->p-1)*R[t][t]/w->Tp[t];
	
	mu[t] = w->W22_0+w->rp/w->Tp[t]*w->p*0.5*pow(Cbar[t],w->p)+w->r2/w->T2[t]*Cbar[t]*Cbar[t];
	mu[t]+= h*(w->Wp2/w->Tp[t]*auxp+w->W22/w->T2[t]*aux2);

	return 0;
}

float compute_Cbar(float*mu, float*Cbar, float**C,float** R, int t, psys *w,float h){
	float Cbar_new;
	float aux2,auxp;
	int m;

	aux2=0.5*R[t][0]*Cbar[0]/w->T2[0];
	auxp=0.5*w->p*(w->p-1)*0.5*R[t][0]*pow(C[t][0],w->p-2)*Cbar[0]/w->Tp[0];
	for(m=1;m<=t-1;m++){
		aux2+=R[t][m]*Cbar[m]/w->T2[m];
		auxp+=w->p*(w->p-1)*0.5*R[t][m]*pow(C[t][m],w->p-2)*Cbar[m]/w->Tp[m];
	}
	aux2+=0.5*R[t][t]*Cbar[t]/w->T2[t];
	auxp+=0.5*w->p*(w->p-1)*0.5*R[t][t]*pow(C[t][t],w->p-2)*Cbar[t]/w->Tp[t];

	Cbar_new = Cbar[t]+h*(-mu[t]*Cbar[t]+w->r2/w->T2[t]*Cbar[t]+w->rp/w->Tp[t]*w->p*0.5*pow(Cbar[t],w->p-1));
	Cbar_new+= h*h*(w->W22/w->T2[t]*aux2+w->Wp2/w->Tp[t]*auxp);

	return Cbar_new;
}

float compute_C(int t, int l,float*mu, float*Cbar, float**C,float** R, psys *w,float h){
	int m,n;
	float auxp1,auxp2,aux21,aux22;
	float Cnew;
	float kroneker;

	if(t==l)kroneker=1;
	else kroneker=0;


	aux21 = .5*R[t][0]*C[l][0]/w->T2[0];
	auxp1 = .5*w->p*(w->p-1)*0.5*pow(C[t][0],w->p-2)*R[t][0]*C[l][0]/w->Tp[0];
	// C is filled by the algorithm only when t'<t, in this case we want
	// to swith l and m, so we have to order C.
	for(m=1; m<=l-1;m++){
		aux21 += R[t][m]*C[l][m]/w->T2[m];
		auxp1 += w->p*(w->p-1)*0.5*pow(C[t][m],w->p-2)*R[t][m]*C[l][m]/w->Tp[m];
	}
	for(m=l; m<=t-1;m++){		
		aux21 += R[t][m]*C[m][l]/w->T2[m];
		auxp1 += w->p*(w->p-1)*0.5*pow(C[t][m],w->p-2)*R[t][m]*C[m][l]/w->Tp[m];
	}
	aux21 += .5*R[t][t]*C[t][l]/w->T2[t];
	auxp1 += .5*w->p*(w->p-1)*0.5*pow(C[t][t],w->p-2)*R[t][t]*C[t][l]/w->Tp[t];

	aux22 = .5*C[t][0]*R[l][0]/w->T2[0];
	auxp2 = .5*w->p*0.5*pow(C[t][0],w->p-1)*R[l][0]/w->Tp[0];
	for(n=1;n<=l-1;n++){
		aux22 += C[t][n]*R[l][n]/w->T2[n];
		auxp2 += w->p*0.5*pow(C[t][n],w->p-1)*R[l][n]/w->Tp[n];
	} 
	aux22 += .5*C[t][l]*R[l][l]/w->T2[t];
	auxp2 += .5*w->p*0.5*pow(C[t][l],w->p-1)*R[l][l]/w->Tp[t];

	Cnew = C[t][l]+h*(-mu[t]*C[t][l]+w->r2/w->T2[t]*Cbar[l]*Cbar[t]+w->rp/w->Tp[t]*w->p*0.5*Cbar[l]*pow(Cbar[t],w->p-1));
	Cnew+= h*(w->W22/w->T2[t]*(aux21+aux22)*h+w->Wp2/w->Tp[t]*(auxp1+auxp2)*h);

	return Cnew;
}

float compute_R(int t, int l,float*mu, float*Cbar, float**C,float** R, psys *w,float h){
	int m,n;
	float kroneker;
	float auxp,aux2;
	
	float Rnew;

	if(t==l)kroneker=1;
	else kroneker=0;
	
	aux2 = R[t][l]*R[l][l]/w->T2[0];
	auxp = w->p*(w->p-1)*0.5*pow(C[t][l],w->p-2)*R[t][l]*R[l][l]/w->Tp[0];
	for(m=l+1; m<=t-1;m++){
		aux2 += R[t][m]*R[m][l]/w->T2[m];
		auxp += w->p*(w->p-1)*0.5*pow(C[t][m],w->p-2)*R[t][m]*R[m][l]/w->Tp[m];
	}
	aux2 += R[t][t]*R[t][l]/w->T2[t];
	auxp += w->p*(w->p-1)*0.5*pow(C[t][t],w->p-2)*R[t][t]*R[t][l]/w->Tp[t];
	
	Rnew = R[t][l]+h*(-mu[t]*R[t][l]+h*w->W22/w->T2[t]*aux2+h*w->Wp2/w->Tp[t]*auxp);

	return Rnew;	
}

float print(float ** C,float ** R,float *mu,float *Cbar,int Nmax, float h, FILE* correlation_file, FILE* response_file, FILE* mu_file,FILE* time_file){
	int i,j;
	
	for(i=0; i<=Nmax;i++){
		
		// if(i%10==0){
		if(1){
			for(j=0;j<Nmax-i; j++){
				
				fprintf(correlation_file, "%f ", C[i+j][j]);
				fprintf(response_file, "%f ", R[i+j][j]);
				fprintf(time_file,"%f", h*(float)(i+j));
			}
			
			fprintf(correlation_file,"\n");
			fprintf(response_file,"\n");
			fprintf(time_file,"\n");
			fprintf(mu_file, "%f \t %f \t %f\n",h*(float)i, Cbar[i], mu[i]);
		}
	}
	return 1.1;
}

float print_mu(float ** C,float ** R,float *mu,float *Cbar,int Nmax, float h, FILE* mu_file){
	int i,j;
	
	for(i=0; i<=Nmax;i++){
		fprintf(mu_file, "%f \t %f \t %f\n",h*(float)i, Cbar[i], mu[i]);
	}
	return 1.1;
}

//	The function modifies T2 and Tp according to a given protocol.
//	The current protocol is the one defined in the paper, T2=1 at all times
//	and Tp reduces gradually in time to the Bayes optimal value.
//	If protocol is set to "false", no protocol will be applied, then T2=Tp=1.
void dumping_policy(psys *w,float C, float tau, int Nmax, float h, bool protocol){
	double t; int i;
	if (protocol){
		for(i=0; i<=Nmax-1; i++){
			t = h*(float)i*w->Delta2;

			if(t>1){
				w->T2[i] = 1.;
				w->Tp[i] = 1.+C/w->Deltap*exp(-t/tau);
			}else{
				w->T2[i] = 1.;
				w->Tp[i] = 1.+C/w->Deltap;
			}
		}
	}else{
		for(i=0; i<=Nmax-1; i++){
			w->T2[i] = 1.; w->Tp[i] = 1.;
		}
	}
	return;
}
