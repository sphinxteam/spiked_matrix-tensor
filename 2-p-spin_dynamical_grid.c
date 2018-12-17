#define thisfile "2-p-spin_dynamical_grid.c"
/*****************************************************
Code based on the idea by [Kim&Latz, 2000], written for the work 
[Berthier, L., Biroli, G., Bouchaud, J. P., Kob, W., Miyazaki, K., & Reichman, D. R. (2007)] 
and  adapted for the (2+p)-planted spin model. The code solves numerically the 
Crisanti-Horner-Sommers-Cugliandolo-Kurchan equations, 
refer to [Sarao Mannelli, S., Biroli, G., Cammarota, C., Krzakala, F., Urbani, P., & Zdeborova, L. (2018)] for details.

written:  08/13/05 by Kuni Miyazaki
updated:  10/12/18 by Stefano Sarao Mannelli
******************************************************/

#define Pi  3.141592653589793  
#define Ntmax 1025
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

typedef struct _psys {
	FILE *fout;
	char dir[100];					//	output directory name
	char name[100];					//	output file name
	double p,Delta2,Deltap,Cbar0;	//	system parameters
	double W22,Wp2,r2,rp;			//	auxiliary variables
} psys;

typedef struct _parr {
	double dt,*mu,dmu,*Cbar;
	double **C, **Q, **M2, **N2, **Mp, **Np;
	double **dCh, **dQh, **dM2h, **dN2h, **dMph, **dNph;
	double **dCv, **dQv, **dM2v, **dN2v, **dMpv, **dNpv;
	double *P2, *Pp, *dCbar;
} parr;

typedef struct _pmct {
	int Nt,rpt,Nt2,Nc,Ntexp,itr;
	int cnt;
	double t,fmin,t0,eps;
} pmct;

//	This functions computes f_p'(x)
double fd1(double x, psys *w) {
	return 0.5*(w->p)*pow(x,w->p-1.0)/w->Delta2;
}

//	This function computes f_p''(x)
double fd2(double x, psys *w) {
	return 0.5*(w->p)*(w->p-1.0)*pow(x,w->p-2.0)/w->Delta2;
}

//	This function computes f_p'''(x)
double fd3(double x, psys *w) {
	return 0.5*(w->p)*(w->p-1.0)*(w->p-2.0)*pow(x,w->p-3.0)/w->Delta2;
}

//	This functions computes f_2'(x)
double f21(double x, psys *w) {
	return x/w->Delta2;
}

//	This function computes f_2''(x)
double f22(double x, psys *w) {
	return 1.0/w->Delta2;
}

//	This function computes the spherical constraint.
double mu_t(pmct *z, parr *x, psys *w, int i) {
	double mu=0;	//	set to 1 once mu is close to its limit value.
	int l;

	for(l=1;l<=i-z->Nc;l++) mu+=w->W22*(x->Q[i][l]-x->Q[i][l-1])*(-1.0*(x->M2[i][l+1]+x->N2[i][l+1]*x->C[i][l+1])+8.0*(x->M2[i][l]+x->N2[i][l]*x->C[i][l])+5.0*(x->M2[i][l-1]+x->N2[i][l-1]*x->C[i][l-1]) )/12.0;
	for(l=1;l<=i-z->Nc;l++) mu+=w->Wp2*(x->Q[i][l]-x->Q[i][l-1])*(-1.0*(x->Mp[i][l+1]+x->Np[i][l+1]*x->C[i][l+1])+8.0*(x->Mp[i][l]+x->Np[i][l]*x->C[i][l])+5.0*(x->Mp[i][l-1]+x->Np[i][l-1]*x->C[i][l-1]) )/12.0;
	mu += 1.0;
	mu += x->dmu - (w->W22*x->M2[i][0] + w->Wp2*x->Mp[i][0])*x->C[i][0];
	mu += - w->r2*x->P2[i]*x->Cbar[i] - w->rp*x->Pp[i]*x->Cbar[i];

	return mu;
}

//	This function computes the integral I1C
double I1C(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	sum2 += x->M2[i][m]*x->C[m][j]-x->M2[i][j]*x->C[j][j]; 
	for(l=m+1;l<=i-1;l++) sum2 += x->dM2v[i][l]*(x->C[l][j]-x->C[l-1][j]);
	sum2 += x->dM2v[i][i]*(-x->C[i-1][j]);
	for(l=j+1;l<=m;l++) sum2 -= (x->M2[i][l]-x->M2[i][l-1])*x->dCh[l][j];

	sump += x->Mp[i][m]*x->C[m][j]-x->Mp[i][j]*x->C[j][j]; 
	for(l=m+1;l<=i-1;l++) sump += x->dMpv[i][l]*(x->C[l][j]-x->C[l-1][j]);
	sump += x->dMpv[i][i]*(-x->C[i-1][j]);
	for(l=j+1;l<=m;l++) sump -= (x->Mp[i][l]-x->Mp[i][l-1])*x->dCh[l][j];

	return w->W22*sum2 + w->Wp2*sump;
}

//	This function computes the integral I2C
double I2C(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	for(l=m+1;l<=i;l++) sum2 += x->dN2v[i][l]*(x->Q[i][l]-x->Q[i][l-1])*(x->C[l][j]+x->C[l-1][j]);
	for(l=j+1;l<=m;l++) sum2 += (x->N2[i][l]+x->N2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCh[l][j];


	for(l=m+1;l<=i;l++) sump += x->dNpv[i][l]*(x->Q[i][l]-x->Q[i][l-1])*(x->C[l][j]+x->C[l-1][j]);
	for(l=j+1;l<=m;l++) sump += (x->Np[i][l]+x->Np[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCh[l][j];

	return 0.5*w->W22*sum2 + 0.5*w->Wp2*sump;
}

//	This function computes the integral I3C
double I3C(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	sum2 += x->M2[i][j]*x->Q[j][j]-x->M2[i][0]*x->Q[j][0]; 
	for(l=1;l<=j;l++) sum2 -= (x->M2[i][l]-x->M2[i][l-1])*x->dQv[j][l];

	sump += x->Mp[i][j]*x->Q[j][j]-x->Mp[i][0]*x->Q[j][0]; 
	for(l=1;l<=j;l++) sump -= (x->Mp[i][l]-x->Mp[i][l-1])*x->dQv[j][l];

	return w->W22*sum2 + w->Wp2*sump;
}

//	This function computes the integral I4C
double I4C(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	for(l=1;l<=j;l++) sum2 += (x->N2[i][l]+x->N2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCv[j][l];

	for(l=1;l<=j;l++) sump += (x->Np[i][l]+x->Np[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCv[j][l];

	return 0.5*w->W22*sum2 + 0.5*w->Wp2*sump;
	//return 0.5*sum2 + 0.5*sump;
}

//	This function computes the integral I1Q
double I1Q(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	sum2 += x->M2[i][m]*x->Q[m][j]-x->M2[i][j]*x->Q[j][j]; 
	for(l=m+1;l<=i-1;l++) sum2 += x->dM2v[i][l]*(x->Q[l][j]-x->Q[l-1][j]);
	sum2 += x->dM2v[i][i]*(-x->Q[i-1][j]);
	for(l=j+1;l<=m;l++) sum2 -= (x->M2[i][l]-x->M2[i][l-1])*x->dQh[l][j];

	sump += x->Mp[i][m]*x->Q[m][j]-x->Mp[i][j]*x->Q[j][j]; 
	for(l=m+1;l<=i-1;l++) sump += x->dMpv[i][l]*(x->Q[l][j]-x->Q[l-1][j]);
	sump += x->dMpv[i][i]*(-x->Q[i-1][j]);
	for(l=j+1;l<=m;l++) sump -= (x->Mp[i][l]-x->Mp[i][l-1])*x->dQh[l][j];

	return w->W22*sum2 + w->Wp2*sump;
}

//	This function computes the integral I2Q
double I2Q(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	for(l=m+1;l<=i;l++) sum2 += x->dN2v[i][l]*(x->Q[i][l]-x->Q[i][l-1])*(2.0-x->Q[l][j]-x->Q[l-1][j]);
	for(l=j+1;l<=m;l++) sum2 += (x->N2[i][l]+x->N2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*(1.0-x->dQh[l][j]);

	for(l=m+1;l<=i;l++) sump += x->dNpv[i][l]*(x->Q[i][l]-x->Q[i][l-1])*(2.0-x->Q[l][j]-x->Q[l-1][j]);
	for(l=j+1;l<=m;l++) sump += (x->Np[i][l]+x->Np[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*(1.0-x->dQh[l][j]);

	return 0.5*w->W22*sum2 + 0.5*w->Wp2*sump;
}

//	This function computes the integral I5Cbar
double I5Cbar(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	for(l=1;l<=i-1;l++) sum2 += x->dM2v[i][l]*(x->Cbar[l]-x->Cbar[l-1]);
	sum2 -= x->dM2v[i][i]*x->Cbar[i-1];

	for(l=1;l<=i-1;l++) sump += x->dMpv[i][l]*(x->Cbar[l]-x->Cbar[l-1]);
	sump -= x->dMpv[i][i]*x->Cbar[i-1];

	return w->W22*sum2 + w->Wp2*sump;
}

//	This function computes the integral I6Cbar
double I6Cbar(parr *x, int i, int j, int m, psys *w) {
	int l;
	double sum2 = 0.0, sump =0.0;

	for(l=1;l<=i;l++) sum2 += (x->N2[i][l]+x->N2[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCbar[l];

	for(l=1;l<=i;l++) sump += (x->Np[i][l]+x->Np[i][l-1])*(x->Q[i][l]-x->Q[i][l-1])*x->dCbar[l];

	return 0.5*w->W22*sum2 + 0.5*w->Wp2*sump;
}

//	This function initializes the arrays.
void initialarray(pmct *z, parr *x, psys *w) {
	int i,j;

	for(i=0;i<=z->Nt2;i++) {
		for(j=0;j<=i;j++) {
			x->C[i][j]= 1.0 - (double)(i-j)*x->dt;
			x->Q[i][j]= 0.0;
			x->M2[i][j]= f21(x->C[i][j],w);
			x->N2[i][j]= f22(x->C[i][j],w);
			x->Mp[i][j]= fd1(x->C[i][j],w);
			x->Np[i][j]= fd2(x->C[i][j],w);
		}
		x->mu[i] = 1.0 - w->W22*f21(x->C[i][0],w)*x->C[i][0] - w->Wp2*fd1(x->C[i][0],w)*x->C[i][0];
		x->Cbar[i] = w->Cbar0+((w->Cbar0*w->Cbar0-1.0)*(w->r2*f21(w->Cbar0,w)+w->rp*fd1(w->Cbar0,w))-w->Cbar0)*i*x->dt;
		x->P2[i] = f21(x->Cbar[i],w);
		x->Pp[i] = fd1(x->Cbar[i],w);
	}
	for(i=1;i<=z->Nt2;i++) {
		for(j=0;j<=i-1;j++) {
			x->dCh[i][j]= 0.5*(x->C[i][j]+x->C[i-1][j]);
			x->dQh[i][j]= 0.5*(x->Q[i][j]+x->Q[i-1][j]);
			x->dM2h[i][j]= 0.5*(x->M2[i][j]+x->M2[i-1][j]);
			x->dN2h[i][j]= 0.5*(x->N2[i][j]+x->N2[i-1][j]);
			x->dMph[i][j]= 0.5*(x->Mp[i][j]+x->Mp[i-1][j]);
			x->dNph[i][j]= 0.5*(x->Np[i][j]+x->Np[i-1][j]);
		}
		x->dCbar[i]= 0.5*(x->Cbar[i]+x->Cbar[i-1]);
	}
	for(i=1;i<=z->Nt2;i++) {
		for(j=1;j<=i;j++) {
			x->dCv[i][j]= 0.5*(x->C[i][j]+x->C[i][j-1]);
			x->dQv[i][j]= 0.5*(x->Q[i][j]+x->Q[i][j-1]);
			x->dM2v[i][j]= 0.5*(x->M2[i][j]+x->M2[i][j-1]);
			x->dN2v[i][j]= 0.5*(x->N2[i][j]+x->N2[i][j-1]);
			x->dMpv[i][j]= 0.5*(x->Mp[i][j]+x->Mp[i][j-1]);
			x->dNpv[i][j]= 0.5*(x->Np[i][j]+x->Np[i][j-1]);
		}
	}
}

//	This function updates quantities in self-consistent loop and returns the error.
double SC(double *gC, double *gQ, double *Cbar, double D, double mu, pmct *z, parr *x, psys *w, int i, int j) {
	int m;
	double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q,i5Cbar,i6Cbar, gCbar=0.0;

	m = (i+j)>>1;
	i1C = I1C(x,i,j,m,w);
	i2C = I2C(x,i,j,m,w);
	i3C = I3C(x,i,j,m,w);
	i4C = I4C(x,i,j,m,w);
	i1Q = I1Q(x,i,j,m,w);
	i2Q = I2Q(x,i,j,m,w);
	i3Q = i3C;
	i4Q = i4C;
	if(j==0){
		i5Cbar = I5Cbar(x,i,j,m,w);
		i6Cbar = I6Cbar(x,i,j,m,w);
  	}

	gC[j] = -0.5/x->dt*x->C[i-2][j]+2.0/x->dt*x->C[i-1][j]-i1C+i2C+i3C+i4C;
	gC[j]+= - w->W22*x->M2[i][0]*x->C[j][0] - w->Wp2*x->Mp[i][0]*x->C[j][0];
	gC[j]+= - (w->r2*x->P2[i] + w->rp*x->Pp[i])*x->Cbar[j];
	gC[j]/= D;
	gC[j]-= x->C[i][j];

	gQ[j] = -1.0 +mu -0.5/x->dt*x->Q[i-2][j] + 2.0/x->dt*x->Q[i-1][j]-i1Q-i2Q-i3Q-i4Q;
	gQ[j]+= w->W22*x->M2[i][0]*x->C[j][0] + w->Wp2*x->Mp[i][0]*x->C[j][0];
	gQ[j]+= + (w->r2*x->P2[i] + w->rp*x->Pp[i])*x->Cbar[j];
	gQ[j]/= D;
	gQ[j]-= x->Q[i][j];

	if(j==0){
		gCbar = -0.5/x->dt*x->Cbar[i-2] + 2.0/x->dt*x->Cbar[i-1];
		gCbar+= - w->r2*x->P2[i] - w->rp*x->Pp[i];
		gCbar+= -i5Cbar + i6Cbar - w->W22*x->M2[i][0]*x->Cbar[0] - w->Wp2*x->Mp[i][0]*x->Cbar[0];
		gCbar/= D;
		x->Cbar[i] = gCbar;
		gCbar = abs(gCbar-Cbar[i-1]);
		return sqrt(gC[j]*gC[j]+gQ[j]*gQ[j]+gCbar*gCbar);
	}

	return sqrt(gC[j]*gC[j]+gQ[j]*gQ[j]);
}

//	This function evaluates all the quantities at a given t.
int step(int i, pmct *z, parr *x, psys *w, int itr) {
	int j,scmax,jmax=0;
	double D,err,emax;
	static double gC[Ntmax],gQ[Ntmax];

// (1) copy the value for the top Nt/4 from the previous column.
	for(j=i-z->Nc;j<=i;j++) {
		x->C[i][j]= x->C[i-1][j-1];
		x->Q[i][j]= x->Q[i-1][j-1];
		x->M2[i][j]= x->M2[i-1][j-1];
		x->N2[i][j]= x->N2[i-1][j-1];
		x->Mp[i][j]= x->Mp[i-1][j-1];
		x->Np[i][j]= x->Np[i-1][j-1];
	}
	x->Cbar[i] = x->Cbar[i-1];
	x->P2[i] = x->P2[i-1];
	x->Pp[i] = x->Pp[i-1];
	for(j=i-z->Nc;j<=i-1;j++) {
		x->dCh[i][j]= x->dCh[i-1][j-1];
		x->dQh[i][j]= x->dQh[i-1][j-1];    
		x->dM2h[i][j]= x->dM2h[i-1][j-1];
		x->dN2h[i][j]= x->dN2h[i-1][j-1];
		x->dMph[i][j]= x->dMph[i-1][j-1];
		x->dNph[i][j]= x->dNph[i-1][j-1];
	}
	for(j=i-z->Nc+1;j<=i;j++) {
		x->dCv[i][j]= x->dCv[i-1][j-1];
		x->dQv[i][j]= x->dQv[i-1][j-1];
		x->dM2v[i][j]= x->dM2v[i-1][j-1];
		x->dN2v[i][j]= x->dN2v[i-1][j-1];
		x->dMpv[i][j]= x->dMpv[i-1][j-1];
		x->dNpv[i][j]= x->dNpv[i-1][j-1];
	}
	x->dCbar[i] = x->dCbar[i-1];

// (2) Prepare test values.
	for(j=0;j<=i-z->Nc-1;j++) {
		x->C[i][j]= x->C[i-1][j];
		x->Q[i][j]= x->Q[i-1][j];
		x->M2[i][j]= f21(x->C[i][j],w);
		x->N2[i][j]= f22(x->C[i][j],w);
		x->Mp[i][j]= fd1(x->C[i][j],w);
		x->Np[i][j]= fd2(x->C[i][j],w);
	}
	for(j=0;j<=i-z->Nc-1;j++) {
		x->dCh[i][j]=(-x->C[i-2][j]+8.0*x->C[i-1][j]+5.0*x->C[i][j])/12.0;
		x->dQh[i][j]=(-x->Q[i-2][j]+8.0*x->Q[i-1][j]+5.0*x->Q[i][j])/12.0;
		x->dM2h[i][j]=(-x->M2[i-2][j]+8.0*x->M2[i-1][j]+5.0*x->M2[i][j])/12.0;
		x->dN2h[i][j]=(-x->N2[i-2][j]+8.0*x->N2[i-1][j]+5.0*x->N2[i][j])/12.0;
		x->dMph[i][j]=(-x->Mp[i-2][j]+8.0*x->Mp[i-1][j]+5.0*x->Mp[i][j])/12.0;
		x->dNph[i][j]=(-x->Np[i-2][j]+8.0*x->Np[i-1][j]+5.0*x->Np[i][j])/12.0;
	}
	for(j=1;j<=i-z->Nc;j++) {
		x->dCv[i][j]= (-x->C[i][j+1]+8.0*x->C[i][j]+5.0*x->C[i][j-1])/12.0;
		x->dQv[i][j]= (-x->Q[i][j+1]+8.0*x->Q[i][j]+5.0*x->Q[i][j-1])/12.0;
		x->dM2v[i][j]= (-x->M2[i][j+1]+8.0*x->M2[i][j]+5.0*x->M2[i][j-1])/12.0;
		x->dN2v[i][j]= (-x->N2[i][j+1]+8.0*x->N2[i][j]+5.0*x->N2[i][j-1])/12.0;
		x->dMpv[i][j]= (-x->Mp[i][j+1]+8.0*x->Mp[i][j]+5.0*x->Mp[i][j-1])/12.0;
		x->dNpv[i][j]= (-x->Np[i][j+1]+8.0*x->Np[i][j]+5.0*x->Np[i][j-1])/12.0;
	}
	x->dCbar[i] = (-x->Cbar[i-2]+8.0*x->Cbar[i-1]+5.0*x->Cbar[i])/12.0;
	x->mu[i] = mu_t(z,x,w,i);
	D = 1.5/x->dt  + x->mu[i] + w->W22*x->dM2v[i][i] + w->Wp2*x->dMpv[i][i];

// (3) Go to the SC loop.
	emax = err =1.0;
	scmax = 0;
	while( err >= z->eps && scmax < z->rpt) {
		for(j=i-z->Nc-1;j>=0;j--) {
			err = SC(gC,gQ,x->Cbar,D,x->mu[i],z,x,w,i,j);
// renew all variable.
			x->C[i][j]+= gC[j];
			if(x->C[i][j] >= z->fmin) x->Q[i][j]+= gQ[j];
			else x->C[i][j] = 0.0;
			x->M2[i][j]= f21(x->C[i][j],w);
			x->N2[i][j]= f22(x->C[i][j],w);
			x->Mp[i][j]= fd1(x->C[i][j],w);
			x->Np[i][j]= fd2(x->C[i][j],w);
			if(j!=i-z->Nc) {
				x->dCh[i][j]=(-x->C[i-2][j]+8.0*x->C[i-1][j]+5.0*x->C[i][j])/12.0;
				x->dQh[i][j]=(-x->Q[i-2][j]+8.0*x->Q[i-1][j]+5.0*x->Q[i][j])/12.0;
				x->dM2h[i][j]=(-x->M2[i-2][j]+8.0*x->M2[i-1][j]+5.0*x->M2[i][j])/12.0;
				x->dN2h[i][j]=(-x->N2[i-2][j]+8.0*x->N2[i-1][j]+5.0*x->N2[i][j])/12.0;
				x->dMph[i][j]=(-x->Mp[i-2][j]+8.0*x->Mp[i-1][j]+5.0*x->Mp[i][j])/12.0;
				x->dNph[i][j]=(-x->Np[i-2][j]+8.0*x->Np[i-1][j]+5.0*x->Np[i][j])/12.0;
			}
			if(j!=0) {
				x->dCv[i][j]= (-x->C[i][j+1]+8.0*x->C[i][j]+5.0*x->C[i][j-1])/12.0;
				x->dQv[i][j]= (-x->Q[i][j+1]+8.0*x->Q[i][j]+5.0*x->Q[i][j-1])/12.0;
				x->dM2v[i][j]= (-x->M2[i][j+1]+8.0*x->M2[i][j]+5.0*x->M2[i][j-1])/12.0;
				x->dN2v[i][j]= (-x->N2[i][j+1]+8.0*x->N2[i][j]+5.0*x->N2[i][j-1])/12.0;
				x->dMpv[i][j]= (-x->Mp[i][j+1]+8.0*x->Mp[i][j]+5.0*x->Mp[i][j-1])/12.0;
				x->dNpv[i][j]= (-x->Np[i][j+1]+8.0*x->Np[i][j]+5.0*x->Np[i][j-1])/12.0;
			}
			x->dCbar[i] = (-x->Cbar[i-2]+8.0*x->Cbar[i-1]+5.0*x->Cbar[i])/12.0;
			x->P2[i] = f21(x->Cbar[i],w);
			x->Pp[i] = fd1(x->Cbar[i],w);
			x->mu[i] = mu_t(z,x,w,i);
			D = 1.5/x->dt  + x->mu[i] + w->W22*x->dM2v[i][i] + w->Wp2*x->dMpv[i][i];
		}
		scmax++;
	}
	return scmax;
}

//	This function halves the values retained and doubles the grid.
void contract(psys *w, pmct *z, parr *x, double *dt, double *dmu) {
	int i,j;
	double Dl=0;

	i=z->Nt;
	for(j=z->Nt-(z->Nc)*2+1;j<=z->Nt-z->Nc;j++) {
		Dl =w->W22*(x->Q[i][j]-x->Q[i][j-1])*(-1.0*(x->M2[i][j+1]+x->N2[i][j+1]*x->C[i][j+1])+8.0*(x->M2[i][j]+x->N2[i][j]*x->C[i][j])+5.0*(x->M2[i][j-1]+x->N2[i][j-1]*x->C[i][j-1]) )/12.0;
		Dl+=w->Wp2*(x->Q[i][j]-x->Q[i][j-1])*(-1.0*(x->Mp[i][j+1]+x->Np[i][j+1]*x->C[i][j+1])+8.0*(x->Mp[i][j]+x->Np[i][j]*x->C[i][j])+5.0*(x->Mp[i][j-1]+x->Np[i][j-1]*x->C[i][j-1]) )/12.0;
		(*dmu) += Dl;
	}
	for(i=1;i<=z->Nt2;i++) {
		for(j=0;j<=i-1;j++) {
			x->dCh[i][j]= 0.5*(x->dCh[2*i][2*j]+x->dCh[2*i-1][2*j]);
			x->dQh[i][j]= 0.5*(x->dQh[2*i][2*j]+x->dQh[2*i-1][2*j]);
			x->dM2h[i][j]= 0.5*(x->dM2h[2*i][2*j]+x->dM2h[2*i-1][2*j]);
			x->dN2h[i][j]= 0.5*(x->dN2h[2*i][2*j]+x->dN2h[2*i-1][2*j]);
			x->dMph[i][j]= 0.5*(x->dMph[2*i][2*j]+x->dMph[2*i-1][2*j]);
			x->dNph[i][j]= 0.5*(x->dNph[2*i][2*j]+x->dNph[2*i-1][2*j]);

		}
		for(j=1;j<=i;j++) {
			x->dCv[i][j]= 0.5*(x->dCv[2*i][2*j]+x->dCv[2*i][2*j-1]);
			x->dQv[i][j]= 0.5*(x->dQv[2*i][2*j]+x->dQv[2*i][2*j-1]);
			x->dM2v[i][j]= 0.5*(x->dM2v[2*i][2*j]+x->dM2v[2*i][2*j-1]);
			x->dN2v[i][j]= 0.5*(x->dN2v[2*i][2*j]+x->dN2v[2*i][2*j-1]);
			x->dMpv[i][j]= 0.5*(x->dMpv[2*i][2*j]+x->dMpv[2*i][2*j-1]);
			x->dNpv[i][j]= 0.5*(x->dNpv[2*i][2*j]+x->dNpv[2*i][2*j-1]);

		}
		x->dCbar[i] = 0.5*(x->dCbar[2*i]+x->dCbar[2*i-1]);
	}
	for(i=0;i<=z->Nt2;i++) {
		for(j=0;j<=i;j++) {
			x->C[i][j]= x->C[2*i][2*j];
			x->Q[i][j]= x->Q[2*i][2*j];
			x->M2[i][j]= x->M2[2*i][2*j];
			x->N2[i][j]= x->N2[2*i][2*j];
			x->Mp[i][j]= x->Mp[2*i][2*j];
			x->Np[i][j]= x->Np[2*i][2*j];
		}
		x->Cbar[i]= x->Cbar[2*i];
	}
	(*dt) *= 2.0;
}

//	This function initializes the output files.
int write_open(pmct *z, parr *x, psys *w) {
	int j=1;

	if((w->fout=fopen(w->name, "wt"))==NULL) {
		fprintf(stderr,"Cannot open the output file:%s\n",strerror(errno));
		return errno;
	}

	fprintf(w->fout,"# (%d) scmax[i]\n",j++);
	fprintf(w->fout,"# (%d) itr\n",j++);
	fprintf(w->fout,"# (%d) t\n",j++);
	fprintf(w->fout,"# (%d) mu\n",j++);
	fprintf(w->fout,"# (%d) Cbar[i]\n",j++);

	fprintf(w->fout,"####[Parameters]##########################");
	fprintf(w->fout,"##########################################\n");
	fprintf(w->fout,"# Program name:   '%s'\n",thisfile);
	fprintf(w->fout,"# This file name: '%s'\n",w->name);
	fprintf(w->fout,"####[System-related parameters]############");
	fprintf(w->fout,"##########################################\n");
	fprintf(w->fout,"# w.p      =%.0f  \n",   w->p);
	fprintf(w->fout,"# w.Delta2      =%.5f\n", w->Delta2);
	fprintf(w->fout,"# w.Deltap      =%.5f\n", w->Deltap);
	fprintf(w->fout,"# z.Nt     =2^%d=%d\n",   z->Ntexp,z->Nt);
	fprintf(w->fout,"# z.Nc    =%d    \n",   z->Nc);
	fprintf(w->fout,"# z.itr    =%d    \n",   z->itr);
	fprintf(w->fout,"# z.fmin   =%.1e \n",   z->fmin);
	fprintf(w->fout,"# z.t0     =%.2e \n",   z->t0);
	fprintf(w->fout,"# z.rpt    =%d   \n",   z->rpt);
	fprintf(w->fout,"# z.eps    =%.1e \n",   z->eps);
	fprintf(w->fout,"# z.Nt2    =%d   \n",   z->Nt2);
	fprintf(w->fout,"# z.Nc     =%d   \n",   z->Nc);
	fprintf(w->fout,"##########################################\n");

	return 0;
}

//	This function updates the output file mu.data, and produces new ones.
int write_add(parr *x, pmct *z, psys *w, int *scmax, int ini, int ifi, int itr) {
	FILE *data;
	char name[100];
	int i;

	for(i=ini;i<=ifi;i++) fprintf(w->fout,"%5d\t%5d\t%12.10e\t%12.10e\t%12.10e\n",scmax[i],itr,(double)i*x->dt,x->mu[i],x->Cbar[i]);
	
	sprintf(name,"%s/%d.data",w->dir,itr);

	if ((data=fopen(name, "wt"))==NULL) {
		fprintf(stderr,"Cannot open the output file:%s\n",strerror(errno));
		return errno;
	}

	fprintf(data,"#\t%d\t%12.10e\n",ini,(double)ini*x->dt);
	for(i=ini;i<=ifi;i++) fprintf(data,"%12.10e\t%12.10e\t%12.10e\n",(double)(i-ini)*x->dt/w->Delta2,x->C[i][ini],(1.0-x->C[i][ini]-x->Q[i][ini])/w->Delta2);

	fflush(w->fout);
	fclose(data);

	return 0;
}

//	This function runs the integration.
int mct(pmct *z, parr *x, psys *w) {
	int i,itr=0;
	int ec;			// return error code;
	int *scmax;

	x->mu = (double *)calloc(Ntmax,sizeof(double));
	x->Cbar = (double *)calloc(Ntmax,sizeof(double));
	x->dCbar = (double *)calloc(Ntmax,sizeof(double));
	x->P2 = (double *)calloc(Ntmax,sizeof(double));
	x->Pp = (double *)calloc(Ntmax,sizeof(double));
	x->C  = (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->C[i] = (double *)calloc(Ntmax,sizeof(double));
	x->Q  = (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->Q[i] = (double *)calloc(Ntmax,sizeof(double));
	x->M2  = (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->M2[i] = (double *)calloc(Ntmax,sizeof(double));
	x->N2  = (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->N2[i] = (double *)calloc(Ntmax,sizeof(double));
	x->Mp  = (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->Mp[i] = (double *)calloc(Ntmax,sizeof(double));
	x->Np  = (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->Np[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dCh= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dCh[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dQh= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dQh[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dM2h= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dM2h[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dN2h= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dN2h[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dMph= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dMph[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dNph= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dNph[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dCv= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dCv[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dQv= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dQv[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dM2v= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dM2v[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dN2v= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dN2v[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dMpv= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dMpv[i] = (double *)calloc(Ntmax,sizeof(double));
	x->dNpv= (double **)calloc(Ntmax,sizeof(double *));
	for(i=0;i<Ntmax;i++) x->dNpv[i] = (double *)calloc(Ntmax,sizeof(double));

	scmax = (int *)calloc(Ntmax,sizeof(int));

	for(i=0;i<=z->Nt;i++) scmax[i] = 0;

	ec=write_open(z,x,w);
	if (ec!=0) {
		fprintf(stderr,"%s\n",strerror(ec));
		return ec;
	}

	printf("inititialization:\t");
	// prepare the array btwn 0 <= i,j <= Nt/2
	initialarray(z,x,w);
	write_add(x,z,w,scmax,0,z->Nt2,itr);
	printf("done\n");
	while(itr <= z->itr) {
		printf("loop %d: \t",itr);
		for(i=(z->Nt2)+1;i<=z->Nt;i++) {
			scmax[i]=step(i,z,x,w,itr);
			if (scmax[i]==z->rpt) {
				itr=z->itr;
				break;
			}
		}
		ec=write_add(x,z,w,scmax,(z->Nt2)+1,z->Nt,itr);
		if (ec!=0) {
			fprintf(stderr,"%s\n",strerror(ec));
			itr=z->itr;
		}
		contract(w,z,x,&(x->dt),&(x->dmu));
		printf("done\n");
		itr++;
	}

	fclose(w->fout);

	free(scmax);
	free(x->mu);
	free(x->Cbar);
	free(x->dCbar);
	for(i=0;i<Ntmax;i++) free(x->C[i]);
	free(x->C);
	for(i=0;i<Ntmax;i++) free(x->Q[i]);
	free(x->Q);
	for(i=0;i<Ntmax;i++) free(x->M2[i]);
	free(x->M2);
	for(i=0;i<Ntmax;i++) free(x->N2[i]);
	free(x->N2);
	for(i=0;i<Ntmax;i++) free(x->Mp[i]);
	free(x->Mp);
	for(i=0;i<Ntmax;i++) free(x->Np[i]);
	free(x->Np);
	for(i=0;i<Ntmax;i++) free(x->dCh[i]);
	free(x->dCh);
	for(i=0;i<Ntmax;i++) free(x->dQh[i]);
	free(x->dQh);
	for(i=0;i<Ntmax;i++) free(x->dM2h[i]);
	free(x->dM2h);
	for(i=0;i<Ntmax;i++) free(x->dN2h[i]);
	free(x->dN2h);
	for(i=0;i<Ntmax;i++) free(x->dMph[i]);
	free(x->dMph);
	for(i=0;i<Ntmax;i++) free(x->dNph[i]);
	free(x->dNph);
	for(i=0;i<Ntmax;i++) free(x->dCv[i]);
	free(x->dCv);
	for(i=0;i<Ntmax;i++) free(x->dQv[i]);
	free(x->dQv);
	for(i=0;i<Ntmax;i++) free(x->dM2v[i]);
	free(x->dM2v);
	for(i=0;i<Ntmax;i++) free(x->dN2v[i]);
	free(x->dN2v);
	for(i=0;i<Ntmax;i++) free(x->dMpv[i]);
	free(x->dMpv);
	for(i=0;i<Ntmax;i++) free(x->dNpv[i]);
	free(x->dNpv);

	return ec;
}

//	The function defines the parameters and run mct. 
int main(int argc, char *argv[]) {
	pmct z;
	parr x;
	psys w;
	char cmd[100];
	int ec;
  
	/*	system parameters	*/
	printf("(2+p)-spin p T Ntexp Nc itr dt eps\n");
	w.p	= 3;
	w.Delta2 = .2;
	w.Deltap = 1.;
	w.Cbar0 = 0.0001;

	w.W22 = 1.;
	w.Wp2 = w.Delta2/w.Deltap*2./w.p;
	w.r2 = -1.0;	
	w.rp = -w.Delta2/w.Deltap*2./w.p; 

	/*	simulation parametres		*/
	z.itr  = 100;							//	number of cycle.
	z.Ntexp	= 8;							//	exponential power of the grid size
	z.Nt   = (int)pow(2,(double)z.Ntexp);	//	grid size (2^Ntexp)
	z.Nt2  = (int)z.Nt/2;					//	half-grid
	z.Nc   = 2;								//	grid elements close to diagonal (small tau)
	z.t0	= 2.5e-7;						//	initial time (it is important to put a small number to have enough information at long time)
	z.rpt  = 3000;							//	maximum number of iterations in self-consistent loop
	z.eps	= 1e-13;						//	precision in the self-consistent loop
	z.fmin = 1e-16;							// 	machine zero
	x.dt   = z.t0/z.Nt2;					//	initial time step
	x.dmu  = 0.0;							//	initial incremental value of the lagrange multiplier mu

	sprintf(w.dir,"2+%.0fplanted-Delta2_%.2f-Deltap_%.2f-Cbar%.2f_Nt%d_Nc%dt0%.2e",w.p,w.Delta2,w.Deltap,w.Cbar0,z.Nt,z.Nc,z.t0);
	sprintf(cmd,"mkdir %s",w.dir);
	ec=system(cmd);
	if (ec!=0) {
		fprintf(stderr,"Cannot create the directory: %s\n",strerror(errno));
		return ec;
	
	}
	sprintf(w.name,"%s/mu.data",w.dir);

	ec=mct(&z,&x,&w);

	return ec;
}
