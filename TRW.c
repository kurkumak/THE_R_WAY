#define thisfile "FP_dyn.c"

//INSTRUCTIONS gcc FP_dyn.c lib/share_memory.c -o FP
/*****************************************************
p-spin model. Algorithm for chi_4. Similar to the aging (Kim&Latz, 2000)
written:  08/13/05
upgraded: 03/30/06
******************************************************
Reference article (appendix C): Spontaneous and induced dynamic correlations in glass formers. II. - L.Berthier, G.Biroli, G.-P.Bouchaud, W.Kob, K.Miyazaki, D.R.Reichman 
modified: 22/03/17
******************************************************/

#define Pi  3.141592653589793
#define Ntmax 1 << 20

#define D1 3./2.
#define D2 -2.
#define D3 1./2.
#define I1 5./12.
#define I2 8./12.
#define I3 -1./12.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#define N_PRINT

/**************************PrInT****************************/

#ifndef N_PRINT
#include "lib/share_memory.h"

void print_vector(double *sh, double **R, int l) {

	int j,k;
	double temp_sh;

	temp_sh = *(sh-2);
	*(sh-2) = 0.;
	for(k=0;k<l;k++) {
		for(j=0;j<l;j++) { *(sh+k*l+j) = R[k][j]; }
	}
	*(sh-2) = temp_sh+1.;
}
#endif

/**********************************************************/

struct pmct{
	int Nt,Ntexp,Nt2,Nc,rpt,den,itr;	//Nt: #(grid elements), Ntexp: log2(Nt), rpt, Nt2: Nt/2, Nsh: #(short-grid elements), den: Nt/Nsh , itr: #(cycles)
	double t,Cmin,t0,eps,alpha;
};

struct parr{							//see appendix C
	double R,dt,dmu;
	double *mu,*E;
	double **C,**Q,**f1,**f2;
	double **dCh,**dQh,**df1h,**df2h;
	double **dCv,**dQv,**df1v,**df2v;
};

struct psys{
	char file[100],dir[100];			//file and directory name
	double p,s_eps,s,T,beta;			//p: (p+s)-spin, s_eps: weight of s-interaction, s: (p+s)-spin, T: temperature, beta: inverse temperature of the initial equilibrium
	double *sh,*sh2;
};

struct t_CQ{
	double time;
	double C;
	double Q;
};

void mct(struct pmct z,struct parr *px,struct psys w);
void initialarray(struct pmct z,struct parr *px,struct psys w);
void contract(struct pmct z,struct parr *px,double *dt,double *dmu);
int step(int i,struct pmct z,struct parr *px,struct psys w);
/**************************SeLf-CoNsIsTeNcE****************************/
double SC2(double *gC,double *gQ,double D,double mu,struct pmct z,struct parr x,struct psys w,int i,int j);
double I1C(struct parr x,int i,int j,int m);
double I2C(struct parr x,int i,int j,int m);
double I3C(struct parr x,int i,int j);
double I4C(struct parr x,int i,int j);
double I1Q(struct parr x,int i,int j,int m);
double I2Q(struct parr x,int i,int j,int m);
double power(double x,int p);
double f(double x,struct psys w);
double fd1(double x,struct psys w);
double fd2(double x,struct psys w);
double mu_t(struct pmct z,struct parr x,struct psys w,int i);
double E_t(struct pmct z,struct parr x,struct psys w,int i);
/**************************ArRay****************************/
void array_initialization(struct pmct z,struct parr *px);
/**************************PaRaMeTeRs****************************/
void parameters_initialization(struct pmct *pz,struct parr *px, struct psys *pw, char *argv[]);
/**************************ScReEn****************************/
void save_config(struct pmct z,struct parr x,struct psys w);
void open_config(struct pmct *pz,struct parr *px,struct psys *pw, char *dir);
/**************************OuTpUt****************************/
void write_parameters(struct pmct z,struct parr x,struct psys w);
void write_C_0(struct parr x,struct psys w,int ini,int ifi);
void write_C(struct pmct z,struct parr x,struct psys w,int j);



int main(int argc, char *argv[]){

	#ifdef N_PRINT
		printf("\nWithout printing!\nTo change see #define N_PRINT\n");
	#endif

	if(argc!=5) { printf("Input 4 parameters! $./FP <B> <Be> <eps_s> <log2(Nt)> \n"); exit(0); }

	struct pmct z;
	struct parr x;
	struct psys w;

	parameters_initialization(&z,&x,&w,argv);

	//open_config(&z,&x,&w,w.dir);

	write_parameters(z,x,w);

	/* shared memory */
	#ifndef N_PRINT
		int KEY = z.Ntexp*12;							//  for the printing procedure
		w.sh = start_double_shared_memory(KEY,z.Nt*z.Nt);			//PRINT_INITIALIZE
		w.sh2 = start_double_shared_memory(KEY+1,z.Nt*z.Nt);		//PRINT_INITIALIZE
	#endif

	printf("\nDIRECTORY_OUTPUT: %s\n", w.dir);

	/*run*/
	printf("\n-----------------------------------------------------START-------------------------------------------------------\n");
	printf("------SYSTEM: (%d+eps*%d)-spin glass, eps = %2.3f------------------------------------------------------------------\n",(int)w.p,(int)w.s,w.s_eps);
	printf("------PARAMETERS: T = %1.3f, T' = %1.3f, grid dimension = %d (Nsh = %d), initial time window = %.2e-----\n",w.T,1./w.beta,z.Nt,z.Nc,z.t0);
	printf("-----------------------------------------------------------------------------------------------------------------\n\n");
	mct(z,&x,w);
	printf("\n-------------------------------------------------------END-------------------------------------------------------\n");

	return(0);
}



void mct(struct pmct z,struct parr *px,struct psys w){

/*------------------------------array initialization------------------------------*/
	int i,j,itr;
	int *scmaxx;	//variables #iteration (itr) reached in the self consitence procedure
	int scMAX,i_scMAX;
	
	array_initialization(z,px);

	scmaxx = (int *) calloc ((z.Nt+1),sizeof(int));
	//for(j=0;j<=z.Nt;j++) scmaxx[j] = 0;
/*------------------------------end array initialization------------------------------*/

	itr=0;

	initialarray(z,px,w); 	// prepare the array btwn 0 <= i,j <= Nt/2

	write_C_0(*px,w,0,z.Nt2);
	write_C(z,*px,w,1);

	#ifndef N_PRINT
		print_vector(w.sh,px->C,z.Nt);		//PRINT
		print_vector(w.sh2,px->Q,z.Nt);		//PRINT
	#endif

	while(itr <= z.itr){

		if(scMAX>10) { save_config(z,*px,w); }

		scMAX = 0;
		for(i=z.Nt2+1;i<=z.Nt;i++){
/*------------------------------------------------------------*/
			scmaxx[i]=step(i,z,px,w);	// propagate the solution the array btwn Nt/2 <= i,j <= Nt
			if(scmaxx[i]>scMAX) { scMAX = scmaxx[i]; i_scMAX = i; }
/*------------------------------------------------------------*/
		}

		write_C_0(*px,w,z.Nt2+1,z.Nt);
		write_C(z,*px,w,1);

	#ifndef N_PRINT
		print_vector(w.sh,px->C,z.Nt);		//PRINT
		print_vector(w.sh2,px->Q,z.Nt);		//PRINT
	#endif

/*------------------------------------------------------------*/
		contract(z,px,&(px->dt),&(px->dmu));	// double the size of the system
/*------------------------------------------------------------*/
		printf(" %dth/%d cycle - self_iter %d/%d: %d - window_time: %.2e \n",itr,z.itr,i_scMAX,z.Nt,scMAX,px->dt*z.Nt);
		itr++;
	}

	for(j=2;j<z.Nt;j++) {
		write_C(z,*px,w,j);
	}
}

void initialarray(struct pmct z,struct parr *px,struct psys w){

/*
	px->C[0][0]=1;  px->Q[0][0]=0;
	mu[0]=T+Bp*f1(1.);
	for(i=0;i<=z.Nt2;i++){
		if(i%100==0)printf("%d\n",i);
		px->Q[i][i]=0; px->C[i][i]=1;
		for(j=0;j<=i;j++){
			px->Q[i+1][j]=px->Q[i][j]+h*(-mu[i]*R[i][j]+If2RR(i,j));
			px->Q[j][i+1]=px->Q[i+1][j];
			px->C[i+1][j]=px->C[i][j]+h*(-mu[i]*C[i][j]+If2RC(i,j)+If1R(i,j)+Bp*df1(px->C[i][0])*px->C[j][0]);
			px->C[j][i+1]=px->C[i+1][j];
	}
	//px->Q[i+1][i]=1;
	//px->Q[i][i+1]=px->Q[i+1][i];
	//px->C[i+1][i+1]=1;
	px->mu[i+1]=If1R(i+1,i+1)+If2RC(i+1,i+1)+T+Bp*df1(px->C[i+1][0])*px->C[i+1][0]; // Questa e' la prescrizione migliore per mu

	double ene=-Bp*f(px->C[i][0])-If1R(i,i);
	fprintf(f11,"%f %f %f %f\n",i*h,mu[i],C[i][0],ene);
	}
*/

	int i,j;

	for(i=0;i<=z.Nt2;i++){
		px->mu[i] = 0.0;
		for(j=0;j<=i;j++){
			px->C[i][j]= 1.0 - (double)(i-j)*px->dt*w.T;	// very short time expansion --> time homogeneous initial C
			px->Q[i][j]= 0.0;							// very short time expansion --> FDT initially respected (Q=0)
			px->f1[i][j]= fd1(px->C[i][j],w);
			px->f2[i][j]= fd2(px->C[i][j],w);
		}
	}
	for(i=1;i<=z.Nt2;i++){
		for(j=0;j<i;j++){
			px->dCh[i][j]= 0.5*(px->C[i-1][j]+px->C[i][j]);
			px->dQh[i][j]= 0.5*(px->Q[i-1][j]+px->Q[i][j]);
			px->df1h[i][j]= 0.5*(px->f1[i-1][j]+px->f1[i][j]);
			px->df2h[i][j]= 0.5*(px->f2[i-1][j]+px->f2[i][j]);
			px->dCv[i][j]= 0.5*(px->C[i][j+1]+px->C[i][j]);
			px->dQv[i][j]= 0.5*(px->Q[i][j+1]+px->Q[i][j]);
			px->df1v[i][j]= 0.5*(px->f1[i][j+1]+px->f1[i][j]);
			px->df2v[i][j]= 0.5*(px->f2[i][j+1]+px->f2[i][j]);
		}
	}

}

int step(int i,struct pmct z,struct parr *px,struct psys w){
	int j,scmax;
	double D,err,err_temp;
	double *gC,*gQ;
	gC = (double *) calloc ((z.Nt),sizeof(double));
	gQ = (double *) calloc ((z.Nt),sizeof(double));

	// (1) copy the value for the top Nsh from the previous column
	px->mu[i] = (D1+1.)*px->mu[i-1]+D2*px->mu[i-2]+D3*px->mu[i-3];

	for(j=3;j<=i;j++){
		px->C[i][j]= (D1+1.)*px->C[i-1][j-1]+D2*px->C[i-2][j-2]+D3*px->C[i-3][j-3];
		px->Q[i][j]= (D1+1.)*px->Q[i-1][j-1]+D2*px->Q[i-2][j-2]+D3*px->Q[i-3][j-3];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][j-1]+D2*px->f1[i-2][j-2]+D3*px->f1[i-3][j-3];
		px->f2[i][j]= (D1+1.)*px->f2[i-1][j-1]+D2*px->f2[i-2][j-2]+D3*px->f2[i-3][j-3];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][j-1]+D2*px->dCh[i-2][j-2]+D3*px->dCh[i-3][j-3];
		px->dQh[i][j]= (D1+1.)*px->dQh[i-1][j-1]+D2*px->dQh[i-2][j-2]+D3*px->dQh[i-3][j-3];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][j-1]+D2*px->df1h[i-2][j-2]+D3*px->df1h[i-3][j-3];
		px->df2h[i][j]= (D1+1.)*px->df2h[i-1][j-1]+D2*px->df2h[i-2][j-2]+D3*px->df2h[i-3][j-3];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][j-1]+D2*px->dCv[i-2][j-2]+D3*px->dCv[i-3][j-3];
		px->dQv[i][j]= (D1+1.)*px->dQv[i-1][j-1]+D2*px->dQv[i-2][j-2]+D3*px->dQv[i-3][j-3];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][j-1]+D2*px->df1v[i-2][j-2]+D3*px->df1v[i-3][j-3];
		px->df2v[i][j]= (D1+1.)*px->df2v[i-1][j-1]+D2*px->df2v[i-2][j-2]+D3*px->df2v[i-3][j-3];
	}

	j=2;
		px->C[i][j]= 2.*px->C[i-1][j-1]-px->C[i-2][j-2];
		px->Q[i][j]= 2.*px->Q[i-1][j-1]-px->Q[i-2][j-2];
		px->f1[i][j]= 2.*px->f1[i-1][j-1]-px->f1[i-2][j-2];
		px->f2[i][j]= 2.*px->f2[i-1][j-1]-px->f2[i-2][j-2];
		px->dCh[i][j]= 2.*px->dCh[i-1][j-1]-px->dCh[i-2][j-2];
		px->dQh[i][j]= 2.*px->dQh[i-1][j-1]-px->dQh[i-2][j-2];
		px->df1h[i][j]= 2.*px->df1h[i-1][j-1]-px->df1h[i-2][j-2];
		px->df2h[i][j]= 2.*px->df2h[i-1][j-1]-px->df2h[i-2][j-2];
		px->dCv[i][j]= 2.*px->dCv[i-1][j-1]-px->dCv[i-2][j-2];
		px->dQv[i][j]= 2.*px->dQv[i-1][j-1]-px->dQv[i-2][j-2];
		px->df1v[i][j]= 2.*px->df1v[i-1][j-1]-px->df1v[i-2][j-2];
		px->df2v[i][j]= 2.*px->df2v[i-1][j-1]-px->df2v[i-2][j-2];
	j=1;
		px->C[i][j]= px->C[i-1][j-1];
		px->Q[i][j]= px->Q[i-1][j-1];
		px->f1[i][j]= px->f1[i-1][j-1];
		px->f2[i][j]= px->f2[i-1][j-1];
		px->dCh[i][j]= px->dCh[i-1][j-1];
		px->dQh[i][j]= px->dQh[i-1][j-1];
		px->df1h[i][j]= px->df1h[i-1][j-1];
		px->df2h[i][j]= px->df2h[i-1][j-1];
		px->dCv[i][j]= px->dCv[i-1][j-1];
		px->dQv[i][j]= px->dQv[i-1][j-1];
		px->df1v[i][j]= px->df1v[i-1][j-1];
		px->df2v[i][j]= px->df2v[i-1][j-1];
	j=0;
		px->C[i][j]= 2.*px->C[i-1][0]-px->C[i-2][0];
		px->Q[i][j]= 2.*px->Q[i-1][0]-px->Q[i-2][0];
		px->f1[i][j]= 2.*px->f1[i-1][0]-px->f1[i-2][0];
		px->f2[i][j]= 2.*px->f2[i-1][0]-px->f2[i-2][0];
		px->dCh[i][j]= 2.*px->dCh[i-1][0]-px->dCh[i-2][0];
		px->dQh[i][j]= 2.*px->dQh[i-1][0]-px->dQh[i-2][0];
		px->df1h[i][j]= 2.*px->df1h[i-1][0]-px->df1h[i-2][0];
		px->df2h[i][j]= 2.*px->df2h[i-1][0]-px->df2h[i-2][0];
		px->dCv[i][j]= 2.*px->dCv[i-1][0]-px->dCv[i-2][0];
		px->dQv[i][j]= 2.*px->dQv[i-1][0]-px->dQv[i-2][0];
		px->df1v[i][j]= 2.*px->df1v[i-1][0]-px->df1v[i-2][0];
		px->df2v[i][j]= 2.*px->df2v[i-1][0]-px->df2v[i-2][0];


  // (2) Prepare test values
	/*for(j=0;j<=i-z.Nc-1;j++){
		px->C[i][j]= px->C[i-1][j];
		px->Q[i][j]= px->Q[i-1][j];
		px->f1[i][j]= fd1(px->C[i][j],w);
		px->f2[i][j]= fd2(px->C[i][j],w);
	}
	for(j=0;j<=i-z.Nc-1;j++){
		px->dCh[i][j] = (-1.*px->C[i-2][j]+8.*px->C[i-1][j]+5.*px->C[i][j])/12.;
		px->dQh[i][j] = (-1.*px->Q[i-2][j]+8.*px->Q[i-1][j]+5.*px->Q[i][j])/12.;
		px->df1h[i][j]=(-1.*px->f1[i-2][j]+8.*px->f1[i-1][j]+5.*px->f1[i][j])/12.;
		px->df2h[i][j]=(-1.*px->f2[i-2][j]+8.*px->f2[i-1][j]+5.*px->f2[i][j])/12.;
		px->dCv[i][j] = (-1.*px->C[i][j+2]+8.*px->C[i][j+1]+5.*px->C[i][j])/12.;
		px->dQv[i][j] = (-1.*px->Q[i][j+2]+8.*px->Q[i][j+1]+5.*px->Q[i][j])/12.;
		px->df1v[i][j]=(-1.*px->f1[i][j+2]+8.*px->f1[i][j+1]+5.*px->f1[i][j])/12.;
		px->df2v[i][j]=(-1.*px->f2[i][j+2]+8.*px->f2[i][j+1]+5.*px->f2[i][j])/12.;
	}*/

	// (3) Go to the SC (self-consistence) loop
	scmax = 0;
	err =1.0; 

	j=1;

	while( err >= z.eps && scmax < z.rpt){
	err =0.0;

	//****** PART ----> j<i-1

	while(j>=0){

		D = 1.5/px->dt  + px->mu[i] + px->df1v[i][i-1]/w.T;

		err_temp = SC2(gC,gQ,D,px->mu[i],z,*px,w,i,j);
		if (err_temp>err) { err=err_temp; }
		// renew all variable
		px->C[i][j]+= gC[j]*z.alpha;
		px->Q[i][j]+= gQ[j]*z.alpha;
		//if(x->C[i][j] >= z.Cmin) x->Q[i][j]+= gQ[j];
		//else { x->C[i][j] = 0.0; }
		px->f1[i][j] = fd1(px->C[i][j],w);
		px->f2[i][j] = fd2(px->C[i][j],w);

		px->dCh[i][j]=(-1.*px->C[i-2][j]+8.*px->C[i-1][j]+5.*px->C[i][j])/12.;
		px->dQh[i][j]=(-1.*px->Q[i-2][j]+8.*px->Q[i-1][j]+5.*px->Q[i][j])/12.;
		px->df1h[i][j]=(-1.*px->f1[i-2][j]+8.*px->f1[i-1][j]+5.*px->f1[i][j])/12.;
		px->df2h[i][j]=(-1.*px->f2[i-2][j]+8.*px->f2[i-1][j]+5.*px->f2[i][j])/12.;
		px->dCv[i][j]= (-1.*px->C[i][j+2]+8.*px->C[i][j+1]+5.*px->C[i][j])/12.;
		px->dQv[i][j]= (-1.*px->Q[i][j+2]+8.*px->Q[i][j+1]+5.*px->Q[i][j])/12.;
		px->df1v[i][j]= (-1.*px->f1[i][j+2]+8.*px->f1[i][j+1]+5.*px->f1[i][j])/12.;
		px->df2v[i][j]= (-1.*px->f2[i][j+2]+8.*px->f2[i][j+1]+5.*px->f2[i][j])/12.;

		px->mu[i] = mu_t(z,*px,w,i);

		j--;
	}

	scmax++;
	if(scmax%10==0) { printf("\r%d/%d %d\r",j,i,scmax); fflush(stdout); }
	j=i-z.Nc;
		
	}

	px->E[i] = E_t(z,*px,w,i);

	free(gC);
	free(gQ);
	return scmax;
}

/**************************SeLf-CoNsIsTeNcE****************************/

double SC2(double *gC, double *gQ, double D, double mu, struct pmct z, struct parr x, struct psys w,int i, int j){
  int m;
  double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q;
  m = (int)(0.5*(i+j));
  i1C = I1C(x,i,j,m);
  i2C = I2C(x,i,j,m);
  i3C = I3C(x,i,j);
  i4C = I4C(x,i,j);
  i1Q = I1Q(x,i,j,m);
  i2Q = I2Q(x,i,j,m);
  i3Q = i3C;
  i4Q = i4C;

  gC[j] = -D3/x.dt*x.C[i-2][j]-D2/x.dt*x.C[i-1][j]+(-i1C+i2C+i3C+i4C)/w.T;
  gC[j]-= (1./w.T-w.beta)*fd1(x.C[i][0],w)*x.C[j][0];
  gC[j]/= D;
  gC[j]-= x.C[i][j];
  gQ[j] = -w.T+mu-D3/x.dt*x.Q[i-2][j]-D2/x.dt*x.Q[i-1][j]+(-i1Q-i2Q-i3Q-i4Q)/w.T;
  gQ[j]+= (1./w.T-w.beta)*fd1(x.C[i][0],w)*x.C[j][0];
  gQ[j]/= D;
  gQ[j]-= x.Q[i][j];
	
  return sqrt(gC[j]*gC[j]+gQ[j]*gQ[j]);
}

double I1C(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][m]*x.C[m][j]-x.f1[i][j]*x.C[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x.df1v[i][l-1]*(x.C[l][j]-x.C[l-1][j]);
  }
  sum += x.df1v[i][i-1]*(-x.C[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dCh[l][j];
  }
  return sum;
}

double I2C(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x.df2v[i][l-1]*(x.Q[i][l]-x.Q[i][l-1])*(x.C[l][j]+x.C[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x.f2[i][l]+x.f2[i][l-1])*(x.Q[i][l]-x.Q[i][l-1])*x.dCh[l][j];
  return 0.5*sum;
}

double I3C(struct parr x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][j]*x.Q[j][j]-x.f1[i][0]*x.Q[j][0];
  for(l=1;l<=j;l++) sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dQv[j][l-1];
  return sum;
}

double I4C(struct parr x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  for(l=1;l<=j;l++) sum += (x.f2[i][l]+x.f2[i][l-1])*(x.Q[i][l]-x.Q[i][l-1])*x.dCv[j][l-1];
  return 0.5*sum;
}

double I1Q(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][m]*x.Q[m][j]-x.f1[i][j]*x.Q[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x.df1v[i][l-1]*(x.Q[l][j]-x.Q[l-1][j]);
  }
  sum += x.df1v[i][i-1]*(-x.Q[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dQh[l][j];
  }
  return sum;
}

double I2Q(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x.df2v[i][l-1]*(x.Q[i][l]-x.Q[i][l-1])*(2.0-x.Q[l][j]-x.Q[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x.f2[i][l]+x.f2[i][l-1])*(x.Q[i][l]-x.Q[i][l-1])*(1.0-x.dQh[l][j]);
  return 0.5*sum;
}

void contract(struct pmct z,struct parr *x,double *dt,double *dmu){
  int i,j;
  double Dl;
  i=z.Nt;
  for(j=z.Nt-z.Nc*2+1;j<=z.Nt-z.Nc;j++){
    Dl=(x->Q[i][j]-x->Q[i][j-1])*(I3*(x->f1[i][j+1]+x->f2[i][j+1]*x->C[i][j+1])
				+I2*(x->f1[i][j  ]+x->f2[i][j  ]*x->C[i][j  ])
				+I1*(x->f1[i][j-1]+x->f2[i][j-1]*x->C[i][j-1]));
    (*dmu) += Dl;
  }
  for(i=1;i<=z.Nt2;i++){
    for(j=0;j<=i-1;j++){
      x->dCh[i][j]= 0.5*(x->dCh[2*i][2*j]+x->dCh[2*i-1][2*j]);
      x->dQh[i][j]= 0.5*(x->dQh[2*i][2*j]+x->dQh[2*i-1][2*j]);
      x->df1h[i][j]= 0.5*(x->df1h[2*i][2*j]+x->df1h[2*i-1][2*j]);
      x->df2h[i][j]= 0.5*(x->df2h[2*i][2*j]+x->df2h[2*i-1][2*j]);
      x->dCv[i][j]= 0.5*(x->dCv[2*i][2*j+1]+x->dCv[2*i][2*j]);
      x->dQv[i][j]= 0.5*(x->dQv[2*i][2*j+1]+x->dQv[2*i][2*j]);
      x->df1v[i][j]= 0.5*(x->df1v[2*i][2*j+1]+x->df1v[2*i][2*j]);
      x->df2v[i][j]= 0.5*(x->df2v[2*i][2*j+1]+x->df2v[2*i][2*j]);
    }
  }
  for(i=0;i<=z.Nt2;i++){
    for(j=0;j<=i;j++){
      x->C[i][j]= x->C[2*i][2*j];
      x->Q[i][j]= x->Q[2*i][2*j];
      x->f1[i][j]= x->f1[2*i][2*j];
      x->f2[i][j]= x->f2[2*i][2*j];
    }
  }
  (*dt) *= 2.0;
}

double power(double x,int p){
	double ff=1.;
	int i;
	for(i=0;i<p;i++) { ff *= x; }
	return ff;
}

/*double If2RR(int i, int j){
	double I=0;
	int k;
		for (k=j+1;k<i;k++){
	I+=fd2(C[i][k])*R[i][k]*R[j][k];
	}
	I*=h;
	return I;
}


double If2RC(int i, int j){
	double I=0;
	int k;
	for (k=0;k<i;k++){
		I+=fd2( C[i][k] )*R[i][k]*C[j][k];
	}
	I*=h;
	return I;
}

double If1R(int i, int j){
	double I=0;
	int k;
	for (k=0;k<j;k++){
		I+=fd1( C[i][k] )*R[j][k];
	}
  	I*=h;
	return I;
}*/

double f(double x, struct psys w){
	return 0.5*power(x,w.p) + w.s_eps*0.5*power(x,w.s);
}

double fd1(double x, struct psys w){
	return 0.5*(w.p)*power(x,w.p-1.0) + w.s_eps*0.5*(w.s)*power(x,w.s-1.0);
}

double fd2(double x, struct psys w){
	return 0.5*(w.p)*(w.p-1.0)*power(x,w.p-2.0) + w.s_eps*0.5*(w.s)*(w.s-1.0)*power(x,w.s-2.0);
}

double mu_t(struct pmct z,struct parr x,struct psys w,int i){
	int l;
	double mu;
	mu = w.T + x.dmu;
	for(l=1;l<=i-z.Nc;l++){
		mu+=(x.Q[i][l]-x.Q[i][l-1])*(I3*(x.f1[i][l+1]+x.f2[i][l+1]*x.C[i][l+1])
			+I2*(x.f1[i][l  ]+x.f2[i][l  ]*x.C[i][l  ])
			+I1*(x.f1[i][l-1]+x.f2[i][l-1]*x.C[i][l-1])
		)/w.T;
	}
	mu -= (1./w.T-w.beta)*fd1(x.C[i][0],w)*x.C[i][0];
	return mu;
}

double E_t(struct pmct z,struct parr x,struct psys w,int i){
	int l;
	double E = 0.;
	for(l=1;l<=i;l++){
		E += (x.Q[i][l]-x.Q[i][l-1])*(I3*x.f1[i][l+1]+I2*x.f1[i][l]+I1*x.f1[i][l-1])/w.T;
	}
	E -= w.beta*f(x.C[i][0],w);
	return E;
}

/**************************ArRay****************************/

void array_initialization(struct pmct z,struct parr *px) {
	int i;

	px->mu = (double *) calloc ((z.Nt+1),sizeof(double));
	px->E = (double *) calloc ((z.Nt+1),sizeof(double));
	px->C  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->Q  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->f1  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->f2  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dCh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dQh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df1h= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df2h= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dCv= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dQv= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df1v= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df2v= (double **) malloc ((z.Nt+1) * sizeof(double*));


	for(i=0;i<(z.Nt+1);i++) {
		px->C[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->Q[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->f1[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->f2[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->dCh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dQh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df1h[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df2h[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dCv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dQv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df1v[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df2v[i]= (double *) calloc ((z.Nt+1),sizeof(double));
	}
}

/**************************PaRaMeTeRs****************************/

void parameters_initialization(struct pmct *pz,struct parr *px, struct psys *pw, char *argv[]) {
	/* system parameters */
	pw->T = 1./atof(argv[1]);			//  T
	pw->beta = atof(argv[2]);			//  beta'
	pw->p     = 3.0;					//  f1(x) = x^p/(2T^2)
	pw->s_eps = atof(argv[3]);		//	weight of s-spin
	pw->s     = 4.0;					//  f2(x) = x^s/(2T^2)

	/* mct parameters */
	pz->Ntexp = atoi(argv[4]);         	//  Nt=2^{Ntexp}  [5..10]
	pz->Nt = 1 << pz->Ntexp;
	pz->Nt2= (int)pz->Nt/2;

	pz->den  = 1 << 5;			// Nt/Nc
	pz->Nc = 1 << 2;
	pz->itr = 50;				// number of cycle.

	pz->t0 = 1.0E-4;		// Initial time window
	pz->Cmin = 1.0E-12;		// Minimal value accepted before setting C identically equal to 0.
	pz->eps= 1.0E-12; 		// Accepted error distance between the solution and the discrete equation ->
	pz->rpt= 1000;			// 		in the iteration procedure with z.rpt the maximum number of iterations

	pz->alpha = 1.;		// Coefficient of self-consistence iteration

	/* array parameters */
	px->R = pw->T*pw->beta;		//Out of equilibrium parameter (R=1 --> equilibrium dynamics)
	px->dt = pz->t0/pz->Nt2;	// Grid time set to Initial Grid time
	px->dmu = 0.;				//IMPORTANT: local upgrade of mu

	/* output files */
	sprintf(pw->file,"0.dat");
	sprintf(pw->dir,"CQ_T%.5f_beta%.5fNt%dt0%.2e",pw->T,pw->beta,pz->Nt,pz->t0);
    mkdir(pw->dir, 0700);

}

/**************************ScReEn****************************/

void save_config(struct pmct z,struct parr x,struct psys w) {
	int i,j;
	FILE *f;
	char fn[100];
	sprintf(fn,"%s/config",w.dir);
	if((f=fopen(fn, "w"))==NULL){
		fprintf(stderr,"save_config: Cannot open a outfile\n");
		exit(1);
	}
	//pmct z
	fprintf(f,"%d %d %d %d %d %d %d\n",z.Nt,z.Ntexp,z.Nt2,z.Nc,z.rpt,z.den,z.itr);
	fprintf(f,"%.5e %.5e %.5e %.5e %.5e\n",z.t,z.Cmin,z.t0,z.eps,z.alpha);
	//psys w
	fprintf(f,"%s %s\n",w.file,w.dir);
	fprintf(f,"%.5e %.5e %.5e %.10e %.10e\n",w.p,w.s_eps,w.s,w.T,w.beta);
	//parr x
	fprintf(f,"%.5e %.5e %.5e\n",x.R,x.dt,x.dmu);
	
	for (i=0; i<z.Nt+1; i++) {
		fprintf(f,"%.5e ",x.mu[i]);
	}
	fprintf(f,"\n");

	for (i=0; i<z.Nt+1; i++) {
		for (j=0; j<z.Nt+1; j++) {
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.C[i][j],x.Q[i][j],x.f1[i][j],x.f2[i][j]);
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.dCh[i][j],x.dQh[i][j],x.df1h[i][j],x.df2h[i][j]);
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.dCv[i][j],x.dQv[i][j],x.df1v[i][j],x.df2v[i][j]);
		}	
	}
	fprintf(f,"\n");
}

void open_config(struct pmct *pz,struct parr *px,struct psys *pw, char *dir) {
	int i,j;
	FILE *f;
	char fn[100];
	sprintf(fn,"%s/config",dir);
	if((f=fopen(fn, "r"))==NULL){
		fprintf(stderr,"open_config: Cannot open a outfile\n");
		exit(1);
	}
	//pmct z
	fscanf(f,"%d %d %d %d %d %d %d\n",&pz->Nt,&pz->Ntexp,&pz->Nt2,&pz->Nc,&pz->rpt,&pz->den,&pz->itr);
	fscanf(f,"%lf %lf %lf %lf %lf\n",&pz->t,&pz->Cmin,&pz->t0,&pz->eps,&pz->alpha);
	//psys w
	fscanf(f,"%s %s\n",pw->file,pw->dir);
	fscanf(f,"%lf %lf %lf %lf %lf\n",&pw->p,&pw->s_eps,&pw->s,&pw->T,&pw->beta);
	//parr x
	fscanf(f,"%lf %lf %lf\n",&px->R,&px->dt,&px->dmu);

	array_initialization(*pz,px);
	
	for (i=0; i<pz->Nt+1; i++) {
		fscanf(f,"%lf ",&px->mu[i]);
	}
	fscanf(f,"\n");

	for (i=0; i<pz->Nt+1; i++) {
		for (j=0; j<pz->Nt+1; j++) {
			fscanf(f,"%lf %lf %lf %lf ",&px->C[i][j],&px->Q[i][j],&px->f1[i][j],&px->f2[i][j]);
			fscanf(f,"%lf %lf %lf %lf ",&px->dCh[i][j],&px->dQh[i][j],&px->df1h[i][j],&px->df2h[i][j]);
			fscanf(f,"%lf %lf %lf %lf ",&px->dCv[i][j],&px->dQv[i][j],&px->df1v[i][j],&px->df2v[i][j]);
		}	
	}
	fscanf(f,"\n");
}

/**************************OuTpUt****************************/

void write_parameters(struct pmct z,struct parr x,struct psys w){
  FILE *fout;
  char fn[100];
  sprintf(fn,"%s/par.txt",w.dir);
  if((fout=fopen(fn, "w"))==NULL){
    fprintf(stderr,"write_parameters: Cannot open a outfile\n");
    exit(1);
  }

  fprintf(fout,"####[Parameters]##########################");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# Program name:   '%s'\n",thisfile);
  fprintf(fout,"# This directory name: '%s'\n",w.dir);
  fprintf(fout,"####[System-related parameters]############");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# w.T      =%.5f\n", w.T);
  fprintf(fout,"# w.beta   =%.5f\n", w.beta);
  fprintf(fout,"# z.Nt     =2^%d=%d\n",   z.Ntexp,z.Nt);
  fprintf(fout,"# z.Nc     =%d   \n",   z.Nc);
  fprintf(fout,"# w.p      =%.0f  \n",   w.p);
  fprintf(fout,"# w.s_eps  =%.2e  \n",   w.s_eps);
  fprintf(fout,"# w.s      =%.0f  \n",   w.s);
  fprintf(fout,"# z.itr    =%d    \n",   z.itr);
  fprintf(fout,"# z.Cmin   =%.1e \n",   z.Cmin);
  fprintf(fout,"# z.t0     =%.2e \n",   z.t0);
  fprintf(fout,"# z.rpt    =%d   \n",   z.rpt);
  fprintf(fout,"# z.eps    =%.1e \n",   z.eps);
  fprintf(fout,"# x.R      =w.T*w.beta=%.2e \n",   x.R);
  fprintf(fout,"####[Time stamp]##########################");
  fprintf(fout,"##########################################\n");
  fclose(fout);
}

void write_C_0(struct parr x, struct psys w, int ini, int ifi){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%s",w.dir,w.file);
	if((fout=fopen(fn, "at"))==NULL){
		fprintf(stderr,"write_C_0: Cannot open a outfile\n");
		exit(1);
	}

	for(i=ini;i<=ifi;i++) {
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",(double)i*x.dt,x.C[i][0],x.Q[i][0],x.mu[i],x.E[i]);
  	}
	fclose(fout);
}

void write_C(struct pmct z, struct parr x, struct psys w, int j){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%.2e.dat",w.dir,j*x.dt);
	if((fout=fopen(fn, "w"))==NULL){
		fprintf(stderr," write_C: Cannot open a outfile\n");
		exit(1);
	}

	for(i=1;i<=z.Nt;i++) {
		fprintf(fout,"%1.5e\t%1.5e\n",(double)(i-j)*x.dt,x.C[i][j]);
  	}
	fclose(fout);
}
