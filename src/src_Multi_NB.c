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

void DifNB(double *zi,dataType *dt1, dataType *dt2, dataType *dt3, dataType *dt4,
	   dataType *dt5, dataType *dt6, dataType *dt7, dataType *dt8, dataType *dt9, 
	   int *k,double *sigma, int *ndt)
{
  char *trans="N";
  int incx,incy,i,K;
  double ONE;
  double *eta1, *eta2, *eta3, *eta4,*eta5, *eta6, *eta7, *eta8, *eta9;
  double dif, *tryz, *newmu1, *newmu2, *newmu3, *newmu4, *newmu5, *newmu6, *newmu7, *newmu8, *newmu9;

  eta1 = dvec(*(dt1->p));
  eta2 = dvec(*(dt2->p));
  eta3 = dvec(*(dt3->p));
  eta4 = dvec(*(dt4->p));
  eta5 = dvec(*(dt5->p));
  eta6 = dvec(*(dt6->p));
  eta7 = dvec(*(dt7->p));
  eta8 = dvec(*(dt8->p));
  eta9 = dvec(*(dt9->p));
  
  incx = 1;
  incy = 1;
  ONE = 1.0;
  K=*k;

  tryz = dvec(*k);

  newmu1 = dvec(*(dt1->p));
  newmu2 = dvec(*(dt2->p));
  newmu3 = dvec(*(dt3->p));
  newmu4 = dvec(*(dt4->p));
  newmu5 = dvec(*(dt5->p));
  newmu6 = dvec(*(dt6->p));
  newmu7 = dvec(*(dt7->p));
  newmu8 = dvec(*(dt8->p));
  newmu9 = dvec(*(dt9->p));
 
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

  if((*ndt) > 2){
    dvcopy(eta3,dt3->beta1,*(dt3->p));
    F77_CALL(dgemv)(trans,dt3->p,k,&ONE,dt3->beta2,dt3->p,tryz,&incx,&ONE,eta3,&incy);
    
    for(i=0; i<*(dt3->p);i++){
      newmu3[i] = exp(eta3[i]);
    }

    dif += logLik(dt3->th,newmu3,dt3->xi,dt3->p) -
      logLik(dt3->th,dt3->mu,dt3->xi,dt3->p);
  }

  if((*ndt) > 3){
    dvcopy(eta4,dt4->beta1,*(dt4->p));
    F77_CALL(dgemv)(trans,dt4->p,k,&ONE,dt4->beta2,dt4->p,tryz,&incx,&ONE,eta4,&incy);
    
    for(i=0; i<*(dt4->p);i++){
      newmu4[i] = exp(eta4[i]);
    }

    dif += logLik(dt4->th,newmu4,dt4->xi,dt4->p) -
      logLik(dt4->th,dt4->mu,dt4->xi,dt4->p);
  }

  if((*ndt) > 4){
    dvcopy(eta5,dt5->beta1,*(dt5->p));
    F77_CALL(dgemv)(trans,dt5->p,k,&ONE,dt5->beta2,dt5->p,tryz,&incx,&ONE,eta5,&incy);
    
    for(i=0; i<*(dt5->p);i++){
      newmu5[i] = exp(eta5[i]);
    }

    dif += logLik(dt5->th,newmu5,dt5->xi,dt5->p) -
      logLik(dt5->th,dt5->mu,dt5->xi,dt5->p);
  }

  if((*ndt) > 5){
    dvcopy(eta6,dt6->beta1,*(dt6->p));
    F77_CALL(dgemv)(trans,dt6->p,k,&ONE,dt6->beta2,dt6->p,tryz,&incx,&ONE,eta6,&incy);
    
    for(i=0; i<*(dt6->p);i++){
      newmu6[i] = exp(eta6[i]);
    }

    dif += logLik(dt6->th,newmu6,dt6->xi,dt6->p) -
      logLik(dt6->th,dt6->mu,dt6->xi,dt6->p);
  }

  if((*ndt) > 6){
    dvcopy(eta7,dt7->beta1,*(dt7->p));
    F77_CALL(dgemv)(trans,dt7->p,k,&ONE,dt7->beta2,dt7->p,tryz,&incx,&ONE,eta7,&incy);
    
    for(i=0; i<*(dt7->p);i++){
      newmu7[i] = exp(eta7[i]);
    }

    dif += logLik(dt7->th,newmu7,dt7->xi,dt7->p) -
      logLik(dt7->th,dt7->mu,dt7->xi,dt7->p);
  }

  if((*ndt) > 7){
    dvcopy(eta8,dt8->beta1,*(dt8->p));
    F77_CALL(dgemv)(trans,dt8->p,k,&ONE,dt8->beta2,dt8->p,tryz,&incx,&ONE,eta8,&incy);
    
    for(i=0; i<*(dt8->p);i++){
      newmu8[i] = exp(eta8[i]);
    }

    dif += logLik(dt8->th,newmu8,dt8->xi,dt8->p) -
      logLik(dt8->th,dt8->mu,dt8->xi,dt8->p);
  }

  if((*ndt) > 8){
    dvcopy(eta9,dt9->beta1,*(dt9->p));
    F77_CALL(dgemv)(trans,dt9->p,k,&ONE,dt9->beta2,dt9->p,tryz,&incx,&ONE,eta9,&incy);
    
    for(i=0; i<*(dt9->p);i++){
      newmu9[i] = exp(eta9[i]);
    }

    dif += logLik(dt9->th,newmu9,dt9->xi,dt9->p) -
      logLik(dt9->th,dt9->mu,dt9->xi,dt9->p);
  }

  dif = dif - 0.5*F77_CALL(ddot)(&K,tryz,&incx,tryz,&incy) + 0.5*F77_CALL(ddot)(&K,zi,&incx,zi,&incy); 

  if((dif > 0) | (runif(0,1) < exp(dif))){
    dvcopy(zi, tryz, *k);
  }

  Free(eta1);
  Free(eta2);
  Free(eta3);
  Free(eta4);
  Free(eta5);
  Free(eta6);
  Free(eta7);
  Free(eta8);
  Free(eta9);

  Free(newmu1);
  Free(newmu2);
  Free(newmu3);
  Free(newmu4);
  Free(newmu5);
  Free(newmu6);
  Free(newmu7);
  Free(newmu8);
  Free(newmu9);

  Free(tryz);
}


void nbMcmc(double *meanz,double *lastz,int *burnin,int *draw,
	    double *a0,double *b0, int *x0,double *mu0, double *th0,int *p0,
	    double *a1,double *b1, int *x1,double *mu1, double *th1,int *p1,
	    double *a2,double *b2, int *x2,double *mu2, double *th2,int *p2,
	    double *a3,double *b3, int *x3,double *mu3, double *th3,int *p3,
	    double *a4,double *b4, int *x4,double *mu4, double *th4,int *p4,
	    double *a5,double *b5, int *x5,double *mu5, double *th5,int *p5,
	    double *a6,double *b6, int *x6,double *mu6, double *th6,int *p6,
	    double *a7,double *b7, int *x7,double *mu7, double *th7,int *p7,
	    double *a8,double *b8, int *x8,double *mu8, double *th8,int *p8,
	    int *k, double *sigma,int *ndt)	    
{

  int i, j, ID;
  double *tempz,*newz;
  newz = dvec(*k);
  tempz = dvec(*k);
  
  dataType *dt[9];

  for(i=0; i<9; i++){
    dt[i] = (dataType *)malloc(sizeof(dataType));
    if(dt[i] == NULL){
      error("Error: cannot allocate memory for dt[]\n");
    }
  }

  fillData(dt[0],a0,b0,x0,mu0,th0,p0);
  fillData(dt[1],a1,b1,x1,mu1,th1,p1);
  fillData(dt[2],a2,b2,x2,mu2,th2,p2);
  fillData(dt[3],a3,b3,x3,mu3,th3,p3);
  fillData(dt[4],a4,b4,x4,mu4,th4,p4);
  fillData(dt[5],a5,b5,x5,mu5,th5,p5);
  fillData(dt[6],a6,b6,x6,mu6,th6,p6);
  fillData(dt[7],a7,b7,x7,mu7,th7,p7);
  fillData(dt[8],a8,b8,x8,mu8,th8,p8);

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
    DifNB(tempz,dt[0],dt[1],dt[2],dt[3],dt[4],dt[5],dt[6],dt[7],dt[8],k,sigma,ndt);
  }

  /* sampling */
  for(j=0; j<(*draw); j++){
    DifNB(tempz,dt[0],dt[1],dt[2],dt[3],dt[4],dt[5],dt[6],dt[7],dt[8],k,sigma,ndt);
    dvadd(newz,tempz,*k);
  }
  PutRNGstate();
  dvscale(newz,*k,1.0/(*draw));
  dvcopy(meanz,newz,*k);
  dvcopy(lastz,tempz,*k);
  
  Free(newz);
  Free(tempz);
}
  
