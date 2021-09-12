#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h> 
#include <R_ext/Lapack.h>

typedef struct{
  int *xi;   /* categorical or count variable, should be a px1 vector  */
  double *beta1;
  double *beta2;
  double *mu; /* for linear regression case */
  double *th; /* for linear regression case */
  int *p;  /* p variable */
}dataType;

double *dvec(int len){
  double *v;
  v = (double *)Calloc(len, double);
  if(v == NULL){
    error("Error: fail to allocate memory space.\n");
  }
  return(v);
}

void dvcopy(double *a, double *b, int row){
  int i;
  for(i=0; i<row; i++){
    a[i] = b[i];
  }
}

/* x = x + y */
void dvadd(double * x, double *y, int size) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = x[i] + y[i];
  }
}

/* x = x*alpha */
void dvscale(double * x, int size, double alpha) {
  int i;
  for(i = 0; i < size; i++) {
    x[i] = alpha*x[i];
  }
}

void fillData(dataType *dt,double *beta1,double *beta2,int *xi,double *mu, double *th,int *p)
{
  dt->xi = xi;
  dt->beta1 = beta1;
  dt->beta2 = beta2;
  dt->mu = mu;
  dt->th = th;  
  dt->p = p;
}

/* double logPost(double *zi,double *alpha,double *beta,int *xi,int *p,int*k) */
double logLik(double *th,double *mu,int *xi, int *p)
{
  int i,j;
  double loglike;
  int count[*p]; 
  
  for(j=0; j<(*p); j++){
    if (xi[j]==0){
      count[j] = j*0+1;
    }else {
      count[j] = j*0;
    }}
  
  loglike = 0;
  for(i=0; i<(*p);i++){
    loglike = loglike + lgamma(th[i] + xi[i]) - lgamma(th[i]) - lgamma(xi[i] + 1) + th[i] * log(th[i])+
      xi[i] * log(mu[i] + count[i])-(th[i] + xi[i]) * log(th[i] + mu[i]);
    //    printf (", %d",count[i]);
    //  printf (", %f",xi[i] * log(mu[i] + count));
    //printf (", %f",xi[i] * log(mu[i] + count[i]));
  }
  return(loglike);
}

/*if tryz is accepted, accept = accept + 1 */
/* void DifNB(double *zi,double *th,double *mu, double *beta1, double *beta2, int *xi, int *p, int *k, double *sigma) */

void loglike_sum(double *zi,dataType *dt1, dataType *dt2, int *k,double *sigma, int *ndt)		 
{
  char *trans="N";
  int incx,incy,i,K;
  double ONE;
  double *eta1, *eta2;
  double dif, *tryz, *newmu1, *newmu2;

  eta1 = dvec(*(dt1->p));
  eta2 = dvec(*(dt2->p));
  
  incx = 1;
  incy = 1;
  ONE = 1.0;
  K=*k;

  tryz = dvec(*k);

  newmu1 = dvec(*(dt1->p));
  newmu2 = dvec(*(dt2->p));
 
  for(i=0; i<*k;i++){
    tryz[i] = zi[i] + rnorm(0, *sigma);
  }

  dif = 0;

  if((*ndt) > 0){
    dvcopy(eta1,dt1->beta1,*(dt1->p));
    F77_CALL(dgemv)(trans,dt1->p,k,&ONE,dt1->beta2,dt1->p,tryz,&incx,&ONE,eta1,&incy);
    
    for(i=0; i<*(dt1->p);i++){
      newmu1[i] = exp(eta1[i]);
    }
    
    dif += logLik(dt1->th,newmu1,dt1->xi,dt1->p) -
      logLik(dt1->th,dt1->mu,dt1->xi,dt1->p);
  }

  if((*ndt) > 1){
    dvcopy(eta2,dt2->beta1,*(dt2->p));
    F77_CALL(dgemv)(trans,dt2->p,k,&ONE,dt2->beta2,dt2->p,tryz,&incx,&ONE,eta2,&incy);
    
    for(i=0; i<*(dt2->p);i++){
      newmu2[i] = exp(eta2[i]);
    }    

    dif += logLik(dt2->th,newmu2,dt2->xi,dt2->p) -
      logLik(dt2->th,dt2->mu,dt2->xi,dt2->p);
  }

  dif = dif - 0.5*F77_CALL(ddot)(&K,tryz,&incx,tryz,&incy) + 0.5*F77_CALL(ddot)(&K,zi,&incx,zi,&incy); 

  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(zi, tryz, *k);
  }

  Free(eta1);
  Free(eta2);

  Free(newmu1);
  Free(newmu2);

  Free(tryz);
}


void gibbs(double *meanz,double *lastz,int *burnin,int *draw,
	   double *a0,double *b0, int *x0,double *mu0, double *th0,int *p0,
	   double *a1,double *b1, int *x1,double *mu1, double *th1,int *p1,
	   int *k, double *sigma,int *ndt)	    
{

  int i, j, ID;
  double *tempz,*newz;
  newz = dvec(*k);
  tempz = dvec(*k);
  
  dataType *dt[2];

  for(i=0; i<2; i++){
    dt[i] = (dataType *)malloc(sizeof(dataType));
    if(dt[i] == NULL){
      error("Error: cannot allocate memory for dt[]\n");
    }
  }

  fillData(dt[0],a0,b0,x0,mu0,th0,p0);
  fillData(dt[1],a1,b1,x1,mu1,th1,p1);

  ID = 0;
  for(j=0; j<(*k); j++){
      tempz[j] = lastz[ID]; 
      ID = ID + 1;
    }

  /* Important to use GetRNGstate */
  GetRNGstate();
  /* MCMC sampling */

  /* burn.in + the first draw*/
  for(i=0; i<(*burnin); i++){
    loglike_sum(tempz, dt[0], dt[1], k, sigma, ndt);
  }

  /* sampling */
  for(j=0; j<(*draw); j++){
    loglike_sum(tempz, dt[0], dt[1], k, sigma, ndt);
    dvadd(newz,tempz,*k);
  }
  
  PutRNGstate();
  dvscale(newz,*k,1.0/(*draw));
  dvcopy(meanz,newz,*k);
  dvcopy(lastz,tempz,*k);
  
  Free(newz);
  Free(tempz);
}
  
