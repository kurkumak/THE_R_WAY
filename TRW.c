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
//#include <math.h>
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
	double p,s_eps,s,T,beta;		//p: (p+s)-spin, s_eps: weight of s-interaction, s: (p+s)-spin, T: temperature, beta: inverse temperature of the initial equilibrium
};

struct parr{							//see appendix C
	double dt,dmu;
	double *mu,*E;
	double **C,**R,**f1,**f2R;
	double **dCh,**dRh,**df1h,**df2Rh;
	double **dCv,**dRv,**df1v,**df2Rv;
};

struct psys{
	char file[100],dir[100];			//file and directory name
	double *sh,*sh2;
};

/**************************GeNeRaL****************************/
void mct(struct pmct z,struct parr *px,struct psys w);
void initialarray(struct pmct z,struct parr *px);
int step(int i,struct pmct z,struct parr *px,struct psys w);
/**************************GrIdMaNiPuLaTiOn****************************/
void contract(struct pmct z,struct parr *px,double *dt,double *dmu);
void Ih(struct parr *px, int i, int j);
void Iv(struct parr *px, int i, int j);
void Mirroring(struct parr *px, int i, int j);
void Extrapolation(struct parr *px, int i, int j);
/**************************SeLf-CoNsIsTeNcE****************************/
double SC2(double *gC,double *gR,double D,double mu,struct pmct z,struct parr x,int i,int j);
double If2RR(struct pmct z,struct parr x,int i,int j);
double If2RC(struct pmct z,struct parr x,int i,int j);
double If1R(struct pmct z,struct parr x,int i,int j);
double grid_If2RR(struct pmct z,struct parr x,int i,int j);
double grid_If2RC(struct pmct z,struct parr x,int i,int j);
double grid_If1R(struct pmct z,struct parr x,int i,int j);
double I1C(struct parr x,int i,int j,int m);
double I2C(struct parr x,int i,int j,int m);
double I3C(struct parr x,int i,int j);
double I4C(struct parr x,int i,int j);
double I1R(struct parr x,int i,int j,int m);
double I2R(struct parr x,int i,int j,int m);
double power(double x,int p);
double f(double x,struct pmct z);
double f1(double x,struct pmct z);
double f2(double x,struct pmct z);
double mu_t(struct pmct z,struct parr x,int i);
double E_t(struct pmct z,struct parr x,int i);
/**************************ArRay****************************/
void array_initialization(struct pmct z,struct parr *px);
/**************************PaRaMeTeRs****************************/
void parameters_initialization(struct pmct *pz,struct parr *px,struct psys *pw,char *argv[]);
/**************************ScReEn****************************/
void save_config(struct pmct z,struct parr x,struct psys w);
void open_config(struct pmct *pz,struct parr *px,struct psys *pw,char *dir);
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
	printf("------SYSTEM: (%d+eps*%d)-spin glass, eps = %2.3f------------------------------------------------------------------\n",(int)z.p,(int)z.s,z.s_eps);
	printf("------PARAMETERS: T = %1.3f, T' = %1.3f, grid dimension = %d (Nsh = %d), initial time window = %.2e-----\n",z.T,1./z.beta,z.Nt,z.Nc,z.t0);
	printf("-----------------------------------------------------------------------------------------------------------------\n\n");
	mct(z,&x,w);
	printf("\n-------------------------------------------------------END-------------------------------------------------------\n");

	return(0);
}


/**************************GeNeRaL****************************/

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

	initialarray(z,px); 	// prepare the array btwn 0 <= i,j <= Nt/2

	write_C_0(*px,w,0,z.Nt2);
	write_C(z,*px,w,1);

	#ifndef N_PRINT
		print_vector(w.sh,px->C,z.Nt);		//PRINT
		print_vector(w.sh2,px->R,z.Nt);		//PRINT
	#endif

	while(itr < z.itr){

		//if(scMAX>10) { save_config(z,*px,w); }

		scMAX = 0;
		for(i=z.Nt2;i<z.Nt;i++){
/*------------------------------------------------------------*/
			scmaxx[i+1]=step(i+1,z,px,w);	// propagate the solution the array btwn Nt/2 <= i,j <= Nt
			if(scmaxx[i]>scMAX) { scMAX = scmaxx[i]; i_scMAX = i; }
/*------------------------------------------------------------*/
		}

		write_C_0(*px,w,z.Nt2+1,z.Nt);
		write_C(z,*px,w,1);

	#ifndef N_PRINT
		print_vector(w.sh,px->C,z.Nt);		//PRINT
		print_vector(w.sh2,px->R,z.Nt);		//PRINT
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

void initialarray(struct pmct z,struct parr *px){

	int i,j;

	px->C[0][0]=1;
	px->R[0][0]=0; 
	px->mu[0]=z.T+z.beta*f1(1.,z);

	for(i=0;i<z.Nt2;i++){

		for(j=i;j>=0;j--){

			px->C[i+1][j] = px->C[i][j]+px->dt*(-px->mu[i]*px->C[i][j]+grid_If2RC(z,*px,i,j)+grid_If1R(z,*px,i,j)+z.beta*f1(px->C[i][0],z)*px->C[j][0]);
			px->R[i+1][j] = px->R[i][j]+px->dt*(-px->mu[i]*px->R[i][j]+grid_If2RR(z,*px,i,j));
			px->f1[i+1][j] = f1(px->C[i+1][j],z);
			px->f2R[i+1][j] = f2(px->C[i+1][j],z)*px->R[i+1][j];

			Ih(px,i+1,j);
			Iv(px,i+1,j);

			Mirroring(px,i+1,j);
		}

		px->R[i+1][i]=1;					//THIS PASSAGE IS FUNDAMENTAL FOR THE STABILITY
		px->f2R[i+1][i] = f2(px->C[i+1][i],z)*px->R[i+1][i];	//
		Ih(px,i+1,i);						//
		Iv(px,i+1,i);						//
		Mirroring(px,i+1,i);					//

		px->mu[i+1]=grid_If1R(z,*px,i+1,i+1)+grid_If2RC(z,*px,i+1,i+1)+z.T+z.beta*f1(px->C[i+1][0],z)*px->C[i+1][0]; // Questa e' la prescrizione migliore per mu
		px->E[i+1]=-z.beta*f(px->C[i+1][0],z)-grid_If1R(z,*px,i+1,i+1);

		px->C[i+1][i+1]=1;
		px->R[i][i]=0;
	}
}

int step(int i,struct pmct z,struct parr *px,struct psys w){
	int j,scmax=0;
	double D,err2=1.0,err_temp2;

	double *gC,*gR;
	gC = (double *) calloc ((z.Nt),sizeof(double));
	gR = (double *) calloc ((z.Nt),sizeof(double));

	// (1) First extrapolation to begin the self-consistence loop
	px->mu[i] = (D1+1.)*px->mu[i-1]+D2*px->mu[i-2]+D3*px->mu[i-3];
	Extrapolation(px,i,j);

	// (2) Go to the SC (self-consistence) loop
	j=1;

	while(err2 >= z.eps*z.eps && scmax < z.rpt){

	//****** PART ----> j<i-1
		while(j>=0){

			D = D1*px->dt + px->mu[i] + px->df1v[i][i-1]/z.T;

			err_temp2 = SC2(gC,gR,D,px->mu[i],z,*px,i,j);
			if (err_temp2>err2) { err2=err_temp2; }

			// renew all variable
			px->C[i][j]+= gC[j]*z.alpha;
			px->R[i][j]+= gR[j]*z.alpha;
			px->f1[i][j] = f1(px->C[i][j],z);
			px->f2R[i][j] = f2(px->C[i][j],z)*px->R[i][j];

			Ih(px,i,j);
			Iv(px,i,j);

			Mirroring(px,i,j);

			j--;
		}

		px->mu[i] = mu_t(z,*px,i);
		scmax++;
		if(scmax%10==0) { printf("\r%d/%d %d\r",j,i,scmax); fflush(stdout); }
		j=i-z.Nc;
	}

	px->E[i] = E_t(z,*px,i);

	free(gC);
	free(gR);

	return scmax;
}

/**************************GrIdMaNiPuLaTiOn****************************/

void contract(struct pmct z,struct parr *x,double *dt,double *dmu){
	int i,j;
	double Dl;
	i=z.Nt;
	for(j=z.Nt-z.Nc*2+1;j<=z.Nt-z.Nc;j++){
		Dl=(x->R[i][j]-x->R[i][j-1])*(I3*(x->f1[i][j+1]+x->f2R[i][j+1]*x->C[i][j+1])
						+I2*(x->f1[i][j  ]+x->f2R[i][j  ]*x->C[i][j  ])
						+I1*(x->f1[i][j-1]+x->f2R[i][j-1]*x->C[i][j-1]));
		(*dmu) += Dl;
	}
	for(i=1;i<=z.Nt2;i++){
		for(j=0;j<=i-1;j++){
			x->dCh[i][j]= 0.5*(x->dCh[2*i][2*j]+x->dCh[2*i-1][2*j]);
			x->dRh[i][j]= 0.5*(x->dRh[2*i][2*j]+x->dRh[2*i-1][2*j]);
			x->df1h[i][j]= 0.5*(x->df1h[2*i][2*j]+x->df1h[2*i-1][2*j]);
			x->df2Rh[i][j]= 0.5*(x->df2Rh[2*i][2*j]+x->df2Rh[2*i-1][2*j]);
			x->dCv[i][j]= 0.5*(x->dCv[2*i][2*j+1]+x->dCv[2*i][2*j]);
			x->dRv[i][j]= 0.5*(x->dRv[2*i][2*j+1]+x->dRv[2*i][2*j]);
			x->df1v[i][j]= 0.5*(x->df1v[2*i][2*j+1]+x->df1v[2*i][2*j]);
			x->df2Rv[i][j]= 0.5*(x->df2Rv[2*i][2*j+1]+x->df2Rv[2*i][2*j]);
		}
	}
	for(i=0;i<=z.Nt2;i++){
		for(j=0;j<=i;j++){
			x->C[i][j]= x->C[2*i][2*j];
			x->R[i][j]= x->R[2*i][2*j];
 			x->f1[i][j]= x->f1[2*i][2*j];
			x->f2R[i][j]= x->f2R[2*i][2*j];
		}
	}
	(*dt) *= 2.0;
}

void Ih(struct parr *px, int i, int j) {
	if(i-2>=j) {
		px->dCh[i-1][j]=I1*px->C[i][j]+I2*px->C[i-1][j]+I3*px->C[i-2][j];
		px->dRh[i-1][j]=I1*px->R[i][j]+I2*px->R[i-1][j]+I3*px->R[i-2][j];
		px->df1h[i-1][j]=I1*px->f1[i][j]+I2*px->f1[i-1][j]+I3*px->f1[i-2][j];
		px->df2Rh[i-1][j]=I1*px->f2R[i][j]+I2*px->f2R[i-1][j]+I3*px->f2R[i-2][j];
	} else if(i-1==j) {
		px->dCh[i-1][j]=0.5*px->C[i][j]+0.5*px->C[i-1][j];
		px->dRh[i-1][j]=0.5*px->R[i][j]+0.5*px->R[i-1][j];
		px->df1h[i-1][j]=0.5*px->f1[i][j]+0.5*px->f1[i-1][j];
		px->df2Rh[i-1][j]=0.5*px->f2R[i][j]+0.5*px->f2R[i-1][j];
	}
}

void Iv(struct parr *px, int i, int j) {
	if(j+2<=i) {
		px->dCv[i][j]= I1*px->C[i][j]+I2*px->C[i][j+1]+I3*px->C[i][j+2];
		px->dRv[i][j]= I1*px->R[i][j]+I2*px->R[i][j+1]+I3*px->R[i][j+2];
		px->df1v[i][j]= I1*px->f1[i][j]+I2*px->f1[i][j+1]+I3*px->f1[i][j+2];
		px->df2Rv[i][j]= I1*px->f2R[i][j]+I2*px->f2R[i][j+1]+I3*px->f2R[i][j+2];
	} else if (j+1==i) {
		px->dCv[i][j]= 0.5*px->C[i][j]+0.5*px->C[i][j+1];
		px->dRv[i][j]= 0.5*px->R[i][j]+0.5*px->R[i][j+1];
		px->df1v[i][j]= 0.5*px->f1[i][j]+0.5*px->f1[i][j+1];
		px->df2Rv[i][j]= 0.5*px->f2R[i][j]+0.5*px->f2R[i][j+1];
	}
}

void Mirroring(struct parr *px, int i, int j) {
	px->C[j][i] = px->C[i][j];
	px->R[j][i] = px->R[i][j];
	px->f1[j][i] = px->f1[i][j];
	px->f2R[j][i] = px->f2R[i][j];

	px->dCv[j][i-1] = px->dCh[i-1][j];
	px->dRv[j][i-1] = px->dRh[i-1][j];
	px->df1v[j][i-1] = px->df1h[i-1][j];
	px->df2Rv[j][i-1] = px->df2Rh[i-1][j];

	px->dCh[j][i] = px->dCv[i][j];
	px->dRh[j][i] = px->dRv[i][j];
	px->df1h[j][i] = px->df1v[i][j];
	px->df2Rh[j][i] = px->df2Rv[i][j];
}


void Extrapolation(struct parr *px, int i, int j) {
	if(j-3>=0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][j-1]+D2*px->C[i-2][j-2]+D3*px->C[i-3][j-3];
		px->R[i][j]= (D1+1.)*px->R[i-1][j-1]+D2*px->R[i-2][j-2]+D3*px->R[i-3][j-3];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][j-1]+D2*px->f1[i-2][j-2]+D3*px->f1[i-3][j-3];
		px->f2R[i][j]= (D1+1.)*px->f2R[i-1][j-1]+D2*px->f2R[i-2][j-2]+D3*px->f2R[i-3][j-3];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][j-1]+D2*px->dCh[i-2][j-2]+D3*px->dCh[i-3][j-3];
		px->dRh[i][j]= (D1+1.)*px->dRh[i-1][j-1]+D2*px->dRh[i-2][j-2]+D3*px->dRh[i-3][j-3];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][j-1]+D2*px->df1h[i-2][j-2]+D3*px->df1h[i-3][j-3];
		px->df2Rh[i][j]= (D1+1.)*px->df2Rh[i-1][j-1]+D2*px->df2Rh[i-2][j-2]+D3*px->df2Rh[i-3][j-3];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][j-1]+D2*px->dCv[i-2][j-2]+D3*px->dCv[i-3][j-3];
		px->dRv[i][j]= (D1+1.)*px->dRv[i-1][j-1]+D2*px->dRv[i-2][j-2]+D3*px->dRv[i-3][j-3];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][j-1]+D2*px->df1v[i-2][j-2]+D3*px->df1v[i-3][j-3];
		px->df2Rv[i][j]= (D1+1.)*px->df2Rv[i-1][j-1]+D2*px->df2Rv[i-2][j-2]+D3*px->df2Rv[i-3][j-3];
	} else if (j-2==0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][2]+D2*px->C[i-2][2]+D3*px->C[i-3][2];
		px->R[i][j]= (D1+1.)*px->R[i-1][2]+D2*px->R[i-2][2]+D3*px->R[i-3][2];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][2]+D2*px->f1[i-2][2]+D3*px->f1[i-3][2];
		px->f2R[i][j]= (D1+1.)*px->f2R[i-1][2]+D2*px->f2R[i-2][2]+D3*px->f2R[i-3][2];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][2]+D2*px->dCh[i-2][2]+D3*px->dCh[i-3][2];
		px->dRh[i][j]= (D1+1.)*px->dRh[i-1][2]+D2*px->dRh[i-2][2]+D3*px->dRh[i-3][2];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][2]+D2*px->df1h[i-2][2]+D3*px->df1h[i-3][2];
		px->df2Rh[i][j]= (D1+1.)*px->df2Rh[i-1][2]+D2*px->df2Rh[i-2][2]+D3*px->df2Rh[i-3][2];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][2]+D2*px->dCv[i-2][2]+D3*px->dCv[i-3][2];
		px->dRv[i][j]= (D1+1.)*px->dRv[i-1][2]+D2*px->dRv[i-2][2]+D3*px->dRv[i-3][2];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][2]+D2*px->df1v[i-2][2]+D3*px->df1v[i-3][2];
		px->df2Rv[i][j]= (D1+1.)*px->df2Rv[i-1][2]+D2*px->df2Rv[i-2][2]+D3*px->df2Rv[i-3][2];
	} else if (j-1==0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][1]+D2*px->C[i-2][1]+D3*px->C[i-3][1];
		px->R[i][j]= (D1+1.)*px->R[i-1][1]+D2*px->R[i-2][1]+D3*px->R[i-3][1];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][1]+D2*px->f1[i-2][1]+D3*px->f1[i-3][1];
		px->f2R[i][j]= (D1+1.)*px->f2R[i-1][1]+D2*px->f2R[i-2][1]+D3*px->f2R[i-3][1];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][1]+D2*px->dCh[i-2][1]+D3*px->dCh[i-3][1];
		px->dRh[i][j]= (D1+1.)*px->dRh[i-1][1]+D2*px->dRh[i-2][1]+D3*px->dRh[i-3][1];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][1]+D2*px->df1h[i-2][1]+D3*px->df1h[i-3][1];
		px->df2Rh[i][j]= (D1+1.)*px->df2Rh[i-1][1]+D2*px->df2Rh[i-2][1]+D3*px->df2Rh[i-3][1];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][1]+D2*px->dCv[i-2][1]+D3*px->dCv[i-3][1];
		px->dRv[i][j]= (D1+1.)*px->dRv[i-1][1]+D2*px->dRv[i-2][1]+D3*px->dRv[i-3][1];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][1]+D2*px->df1v[i-2][1]+D3*px->df1v[i-3][1];
		px->df2Rv[i][j]= (D1+1.)*px->df2Rv[i-1][1]+D2*px->df2Rv[i-2][1]+D3*px->df2Rv[i-3][1];
	}  else if (j==0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][0]+D2*px->C[i-2][0]+D3*px->C[i-3][0];
		px->R[i][j]= (D1+1.)*px->R[i-1][0]+D2*px->R[i-2][0]+D3*px->R[i-3][0];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][0]+D2*px->f1[i-2][0]+D3*px->f1[i-3][0];
		px->f2R[i][j]= (D1+1.)*px->f2R[i-1][0]+D2*px->f2R[i-2][0]+D3*px->f2R[i-3][0];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][0]+D2*px->dCh[i-2][0]+D3*px->dCh[i-3][0];
		px->dRh[i][j]= (D1+1.)*px->dRh[i-1][0]+D2*px->dRh[i-2][0]+D3*px->dRh[i-3][0];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][0]+D2*px->df1h[i-2][0]+D3*px->df1h[i-3][0];
		px->df2Rh[i][j]= (D1+1.)*px->df2Rh[i-1][0]+D2*px->df2Rh[i-2][0]+D3*px->df2Rh[i-3][0];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][0]+D2*px->dCv[i-2][0]+D3*px->dCv[i-3][0];
		px->dRv[i][j]= (D1+1.)*px->dRv[i-1][0]+D2*px->dRv[i-2][0]+D3*px->dRv[i-3][0];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][0]+D2*px->df1v[i-2][0]+D3*px->df1v[i-3][0];
		px->df2Rv[i][j]= (D1+1.)*px->df2Rv[i-1][0]+D2*px->df2Rv[i-2][0]+D3*px->df2Rv[i-3][0];
	}
}

/**************************SeLf-CoNsIsTeNcE****************************/

double SC2(double *gC, double *gR, double D, double mu, struct pmct z, struct parr x, int i, int j){
  int m;
  double i1C,i2C,i3C,i4C,i1R,i2R,i3R,i4R;
  m = (int)(0.5*(i+j));
  i1C = I1C(x,i,j,m);
  i2C = I2C(x,i,j,m);
  i3C = I3C(x,i,j);
  i4C = I4C(x,i,j);
  i1R = I1R(x,i,j,m);
  i2R = I2R(x,i,j,m);
  i3R = i3C;
  i4R = i4C;

  gC[j] = -D3/x.dt*x.C[i-2][j]-D2/x.dt*x.C[i-1][j]+(-i1C+i2C+i3C+i4C)/z.T;
  gC[j]-= (1./z.T-z.beta)*f1(x.C[i][0],z)*x.C[j][0];
  gC[j]/= D;
  gC[j]-= x.C[i][j];
  gR[j] = -z.T+mu-D3/x.dt*x.R[i-2][j]-D2/x.dt*x.R[i-1][j]+(-i1R-i2R-i3R-i4R)/z.T;
  gR[j]+= (1./z.T-z.beta)*f1(x.C[i][0],z)*x.C[j][0];
  gR[j]/= D;
  gR[j]-= x.R[i][j];
	
  return gC[j]*gC[j]+gR[j]*gR[j];
}

double If2RR(struct pmct z,struct parr x,int i,int j){
	double I=0;
	int k;
	//I+=0.5*(f2(x.C[i][j],z)*x.R[i][j]*x.R[j][j]+f2(x.C[i][i],z)*x.R[i][i]*x.R[j][i]);
	I+=0.5*(x.f2R[i][j]*x.R[j][j]+x.f2R[i][i]*x.R[j][i]);
	for(k=j+1;k<i;k++){
		//I+=f2(x.C[i][k],z)*x.R[i][k]*x.R[j][k];
		I+=x.f2R[i][k]*x.R[j][k];
	}
	I*=x.dt;
	return I;
}


double If2RC(struct pmct z,struct parr x,int i,int j){
	double I=0;
	int k;
	//I+=0.5*(f2(x.C[i][0],z)*x.R[i][0]*x.C[j][0]+f2(x.C[i][i],z)*x.R[i][i]*x.C[j][i]);
	I+=0.5*(x.f2R[i][0]*x.C[j][0]+x.f2R[i][i]*x.C[j][i]);
	for(k=1;k<i;k++){
		//I+=f2(x.C[i][k],z)*x.R[i][k]*x.C[j][k];
		I+=x.f2R[i][k]*x.C[j][k];
	}
	I*=x.dt;
	return I;
}

double If1R(struct pmct z,struct parr x,int i,int j){
	double I=0;
	int k;
	//I+=0.5*(f1(x.C[i][0],z)*x.R[j][0]+f1(x.C[i][j],z)*x.R[j][j]);
	I+=0.5*(x.f1[i][0]*x.R[j][0]+x.f1[i][j]*x.R[j][j]);
	for(k=1;k<j;k++){
		//I+=f1(x.C[i][k],z)*x.R[j][k];
		I+=x.f1[i][k]*x.R[j][k];
	}
	I*=x.dt;
	return I;
}

double grid_If2RR(struct pmct z,struct parr x,int i,int j){
	double I=0;
	int k;
	for(k=j;k<i;k++){
		I+=x.df2Rv[i][k]*x.dRv[j][k];
	}
	I*=x.dt;
	return I;
}


double grid_If2RC(struct pmct z,struct parr x,int i,int j){
	double I=0;
	int k;
	for(k=0;k<i;k++){
		I+=x.df2Rv[i][k]*x.dCv[j][k];
	}
	I*=x.dt;
	return I;
}

double grid_If1R(struct pmct z,struct parr x,int i,int j){
	double I=0;
	int k;
	for(k=0;k<j;k++){
		I+=x.df1v[i][k]*x.dRv[j][k];
	}
	I*=x.dt;
	return I;
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
  for(l=m+1;l<=i;l++) sum += x.df2Rv[i][l-1]*(x.R[i][l]-x.R[i][l-1])*(x.C[l][j]+x.C[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x.f2R[i][l]+x.f2R[i][l-1])*(x.R[i][l]-x.R[i][l-1])*x.dCh[l][j];
  return 0.5*sum;
}

double I3C(struct parr x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][j]*x.R[j][j]-x.f1[i][0]*x.R[j][0];
  for(l=1;l<=j;l++) sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dRv[j][l-1];
  return sum;
}

double I4C(struct parr x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  for(l=1;l<=j;l++) sum += (x.f2R[i][l]+x.f2R[i][l-1])*(x.R[i][l]-x.R[i][l-1])*x.dCv[j][l-1];
  return 0.5*sum;
}

double I1R(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][m]*x.R[m][j]-x.f1[i][j]*x.R[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x.df1v[i][l-1]*(x.R[l][j]-x.R[l-1][j]);
  }
  sum += x.df1v[i][i-1]*(-x.R[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dRh[l][j];
  }
  return sum;
}

double I2R(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x.df2Rv[i][l-1]*(x.R[i][l]-x.R[i][l-1])*(2.0-x.R[l][j]-x.R[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x.f2R[i][l]+x.f2R[i][l-1])*(x.R[i][l]-x.R[i][l-1])*(1.0-x.dRh[l][j]);
  return 0.5*sum;
}

double power(double x,int p){
	double ff=1.;
	int i;
	for(i=0;i<p;i++) { ff *= x; }
	return ff;
}

double f(double x, struct pmct z){
	return 0.5*power(x,z.p) + z.s_eps*0.5*power(x,z.s);
}

double f1(double x, struct pmct z){
	return 0.5*(z.p)*power(x,z.p-1.0) + z.s_eps*0.5*(z.s)*power(x,z.s-1.0);
}

double f2(double x, struct pmct z){
	return 0.5*(z.p)*(z.p-1.0)*power(x,z.p-2.0) + z.s_eps*0.5*(z.s)*(z.s-1.0)*power(x,z.s-2.0);
}

double mu_t(struct pmct z,struct parr x,int i){
	int l;
	double mu;
	mu = z.T + x.dmu;
	for(l=1;l<=i-z.Nc;l++){
		mu+=(x.R[i][l]-x.R[i][l-1])*(I3*(x.f1[i][l+1]+x.f2R[i][l+1]*x.C[i][l+1])
			+I2*(x.f1[i][l  ]+x.f2R[i][l  ]*x.C[i][l  ])
			+I1*(x.f1[i][l-1]+x.f2R[i][l-1]*x.C[i][l-1])
		)/z.T;
	}
	mu -= (1./z.T-z.beta)*f1(x.C[i][0],z)*x.C[i][0];
	return mu;
}

double E_t(struct pmct z,struct parr x,int i){
	int l;
	double E = 0.;
	for(l=1;l<=i;l++){
		E += (x.R[i][l]-x.R[i][l-1])*(I3*x.f1[i][l+1]+I2*x.f1[i][l]+I1*x.f1[i][l-1])/z.T;
	}
	E -= z.beta*f(x.C[i][0],z);
	return E;
}

/**************************ArRay****************************/

void array_initialization(struct pmct z,struct parr *px) {
	int i;

	px->mu = (double *) calloc ((z.Nt+1),sizeof(double));
	px->E = (double *) calloc ((z.Nt+1),sizeof(double));
	px->C  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->R  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->f1  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->f2R  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dCh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dRh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df1h= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df2Rh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dCv= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dRv= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df1v= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df2Rv= (double **) malloc ((z.Nt+1) * sizeof(double*));


	for(i=0;i<(z.Nt+1);i++) {
		px->C[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->R[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->f1[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->f2R[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->dCh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dRh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df1h[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df2Rh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dCv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dRv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df1v[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df2Rv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
	}
}

/**************************PaRaMeTeRs****************************/

void parameters_initialization(struct pmct *pz,struct parr *px, struct psys *pw, char *argv[]) {
	/* system parameters */
	pz->T = 1./atof(argv[1]);		//  T
	pz->beta = atof(argv[2]);		//  beta'
	pz->p     = 3.0;			//  f1(x) = x^p/2
	pz->s_eps = atof(argv[3]);		//	weight of s-spin
	pz->s     = 4.0;			//  f2(x) = x^s/2

	/* mct parameters */
	pz->Ntexp = atoi(argv[4]);         	//  Nt=2^{Ntexp}  [5..10]
	pz->Nt = 1 << pz->Ntexp;
	pz->Nt2= (int)pz->Nt/2;

	pz->den  = 1 << 5;			// Nt/Nc
	pz->Nc = 1 << 2;
	pz->itr = 0;				// number of cycle.

	pz->t0 = 1.0E0;			// Initial time window
	pz->Cmin = 1.0E-12;			// Minimal value accepted before setting C identically equal to 0.
	pz->eps= 1.0E-12; 			// Accepted error distance between the solution and the discrete equation ->
	pz->rpt= 1000;				// 		in the iteration procedure with z.rpt the maximum number of iterations

	pz->alpha = 1.;				// Coefficient of self-consistence iteration

	/* array parameters */
	px->dt = pz->t0/pz->Nt2;		// Grid time set to Initial Grid time
	px->dmu = 0.;				//IMPORTANT: local upgrade of mu

	/* output files */
	sprintf(pw->file,"0.dat");
	sprintf(pw->dir,"CR_T%.5f_beta%.5fNt%dt0%.2e",pz->T,pz->beta,pz->Nt,pz->t0);
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
	fprintf(f,"%.5e %.5e %.5e %.10e %.10e\n",z.p,z.s_eps,z.s,z.T,z.beta);
	//psys w
	fprintf(f,"%s %s\n",w.file,w.dir);
	//parr x
	fprintf(f,"%.5e %.5e\n",x.dt,x.dmu);
	
	for (i=0; i<z.Nt+1; i++) {
		fprintf(f,"%.5e ",x.mu[i]);
	}
	fprintf(f,"\n");

	for (i=0; i<z.Nt+1; i++) {
		for (j=0; j<z.Nt+1; j++) {
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.C[i][j],x.R[i][j],x.f1[i][j],x.f2R[i][j]);
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.dCh[i][j],x.dRh[i][j],x.df1h[i][j],x.df2Rh[i][j]);
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.dCv[i][j],x.dRv[i][j],x.df1v[i][j],x.df2Rv[i][j]);
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
	fscanf(f,"%lf %lf %lf %lf %lf\n",&pz->p,&pz->s_eps,&pz->s,&pz->T,&pz->beta);
	//psys w
	fscanf(f,"%s %s\n",pw->file,pw->dir);
	//parr x
	fscanf(f,"%lf %lf\n",&px->dt,&px->dmu);

	array_initialization(*pz,px);
	
	for (i=0; i<pz->Nt+1; i++) {
		fscanf(f,"%lf ",&px->mu[i]);
	}
	fscanf(f,"\n");

	for (i=0; i<pz->Nt+1; i++) {
		for (j=0; j<pz->Nt+1; j++) {
			fscanf(f,"%lf %lf %lf %lf ",&px->C[i][j],&px->R[i][j],&px->f1[i][j],&px->f2R[i][j]);
			fscanf(f,"%lf %lf %lf %lf ",&px->dCh[i][j],&px->dRh[i][j],&px->df1h[i][j],&px->df2Rh[i][j]);
			fscanf(f,"%lf %lf %lf %lf ",&px->dCv[i][j],&px->dRv[i][j],&px->df1v[i][j],&px->df2Rv[i][j]);
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
  fprintf(fout,"# z.T      =%.5f\n", z.T);
  fprintf(fout,"# z.beta   =%.5f\n", z.beta);
  fprintf(fout,"# z.Nt     =2^%d=%d\n",   z.Ntexp,z.Nt);
  fprintf(fout,"# z.Nc     =%d   \n",   z.Nc);
  fprintf(fout,"# z.p      =%.0f  \n",   z.p);
  fprintf(fout,"# z.s_eps  =%.2e  \n",   z.s_eps);
  fprintf(fout,"# z.s      =%.0f  \n",   z.s);
  fprintf(fout,"# z.itr    =%d    \n",   z.itr);
  fprintf(fout,"# z.Cmin   =%.1e \n",   z.Cmin);
  fprintf(fout,"# z.t0     =%.2e \n",   z.t0);
  fprintf(fout,"# z.rpt    =%d   \n",   z.rpt);
  fprintf(fout,"# z.eps    =%.1e \n",   z.eps);
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
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",(double)i*x.dt,x.C[i][0],x.R[i][0],x.mu[i],x.E[i]);
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
