#include <R.h>
#include <Rdefines.h>
#include<stdlib.h>
#include<Rmath.h>

#define NVAR 7  

typedef struct {
  double L, M, R, logage, pcage;
} DATA5;



static int orderage(const void *p1, const void *p2)
{
  DATA5 *m1 = (DATA5 *) p1;
  DATA5 *m2 = (DATA5 *) p2;
  
  if(m1->logage < m2->logage) 
    return -1;
  else if(m1->logage > m2->logage)
    return 1;
  else
    return 0;
}


void findrange(double *x, long dim, double vl, double vu, long *l, long *u)
{
  long res, up, low;
  
  if(vl <= x[0])
    res = 0;
  else if(vl > x[dim-1])
    res = -1;
  else
    {
      up = dim-1;
      low = 0;
      while(up-low > 1)
	{
	  res = (up+low)/2;
	  if(x[res] < vl)
	    low = res;
	  else 
	    up = res;
	}
      res = up;
    }
  *l = res;

  if(vu < x[0])
    res = -1;
  else if(vu >= x[dim-1])
    res = dim-1;
  else
    {
      up = dim-1;
      low = *l;
      while(up-low > 1)
	{
	  res = (up+low)/2;
	  if(x[res] < vu)
	    low = res;
	  else 
	    up = res;
	}
      res = low;
    }
  *u = res;
}



SEXP lik4bin(SEXP data, SEXP star, SEXP sigma, SEXP thr, SEXP var, SEXP power, SEXP restringi, SEXP tsp)
{
  double *Pdata, *Psigma, *Pstar, Pres[22], *Rres, *wstar, *age1, *age2;
  double *Teff, *logg, *z, *M, *R, *Dni, *nimax, *logage, *pcage;
  double Vthr, maxL, maxL1, maxL2, lmult, EXP, rpcage;
  long nrow, ncol, count;

  double sq2pi, chi[NVAR], locsigma[NVAR], chi2, mult, L, mass, 
    radius, lt, ltnlog;;
  double sTeffP, sTeffM, time1, time2;
  SEXP res, dm, sel;
  long i, j, nres, nres1, nres2, DIM, start, n, startT, stopT, up, low;
  int ii, norun, nstar, *Psel, *Pvar, restr;
  DATA5 *d, *d1, *d2, *d3, *d4;
  long lb, ub;
  double t_spread; // max diff. in age

  // cast and pointers 
  PROTECT(data = AS_NUMERIC(data));
  PROTECT(star = AS_NUMERIC(star));
  PROTECT(sigma = AS_NUMERIC(sigma));
  PROTECT(thr = AS_NUMERIC(thr));
  PROTECT(var = AS_INTEGER(var));
  PROTECT(power = AS_NUMERIC(power));
  PROTECT(restringi = AS_INTEGER(restringi));
  PROTECT(tsp = AS_NUMERIC(tsp));

  Pdata = NUMERIC_POINTER(data);
  Pstar = NUMERIC_POINTER(star);
  Psigma = NUMERIC_POINTER(sigma);
  Vthr = NUMERIC_VALUE(thr);
  Pvar = INTEGER_POINTER(var);
  EXP = NUMERIC_VALUE(power);
  restr = NUMERIC_VALUE(restringi);
  t_spread = NUMERIC_VALUE(tsp);

  // sqrt ( 2 * pi )
  sq2pi = 2.506628274631000;

  // dataset dimensions
  nrow = INTEGER(GET_DIM(data))[0];
  ncol = INTEGER(GET_DIM(data))[1];

  // column pointers
  // data are column ordered!
  Teff = Pdata;
  logg = Pdata+nrow;
  z = Pdata+2*nrow;
  Dni = Pdata+3*nrow;
  nimax = Pdata+4*nrow;
  M = Pdata+5*nrow;
  R = Pdata+6*nrow;
  logage = Pdata+7*nrow;
  pcage = Pdata+8*nrow;

  // vector for likelihood computations
  // 1 = include; 0 = exclude
  Psel = (int*)malloc(nrow*sizeof(int));

  for(nstar=0;nstar<2;nstar++) 
    {
      for(j=0;j<nrow;j++)
	Psel[j] = 0;

      wstar = &Pstar[(nstar)*9];

      // sigma scaling for Dni,nimax,M,R (it is a % in input)
      for(n=0;n<NVAR;n++)
	locsigma[n] = Psigma[n+NVAR*nstar];
      for(n=3;n<7;n++)
	locsigma[n] *= wstar[n];
      
      mult = 1;
      for(n=0;n<NVAR;n++)
	if(Pvar[n] == 1)
	  mult *= 1.0/(sq2pi * locsigma[n]);
      lmult = log(mult);

      // allowed Teff interval
      sTeffP = wstar[0] + Vthr*locsigma[0];
      sTeffM = wstar[0] - Vthr*locsigma[0];

      // ricerca righe con Teff minima e massima
      findrange(Teff, nrow, sTeffM, sTeffP, &startT, &stopT);
      if(startT == -1 || stopT == -1)
	{
	  free(Psel);
	  UNPROTECT(8);
	  return(R_NilValue);
	}

      // sel computation
      nres = 0;
      for(j=startT;j<=stopT;j++)
	{
	  for(ii=0;ii<NVAR;ii++)
	    chi[ii] = 0;
	  
	  if(Pvar[0] == 1)
	    chi[0] = (Teff[j] - wstar[0])/locsigma[0];
	  if(Pvar[1] == 1)
	    chi[1] = (logg[j] - wstar[1])/locsigma[1];
	  if(Pvar[2] == 1)
	    chi[2] = (z[j] - wstar[2])/locsigma[2];
	  if(Pvar[3] == 1)
	    chi[3] = (Dni[j] - wstar[3])/locsigma[3];
	  if(Pvar[4] == 1)
	    chi[4] = (nimax[j] - wstar[4])/locsigma[4];
	  if(Pvar[5] == 1)
	    chi[5] = (M[j] - wstar[5])/locsigma[5];
	  if(Pvar[6] == 1)
	    chi[6] = (R[j] - wstar[6])/locsigma[6];

	  norun = 0;
	  for(ii=0;ii<NVAR;ii++)
	    {
	      if(fabs(chi[ii]) >= Vthr)
		{
		  norun = 1;
		  break;
		}
	    }
	  
	  if( norun == 0 ) 
	    {
	      chi2 = 0;
	      for(ii=0;ii<NVAR;ii++)
		chi2 += chi[ii]*chi[ii];
	      if( restr == 1 ) 
		{
		  if(sqrt(chi2) <= 3 )
		    {
		      nres++;
		      Psel[j] = 1;
		    }
		}
	      else 
		{
		  nres++;
		  Psel[j] = 1;
		}
	    }
	}
      
      // no data! return
      if(nres == 0) 
	{
	  free(Psel);
	  UNPROTECT(8);
	  return(R_NilValue);
	}
      // init output matrix
      DIM = nres;
      if(nstar == 0)
	{
	  d1 = (DATA5 *)calloc(DIM+1, sizeof(DATA5));
	  d = d1;
	}
      else
	{
	  d2 = (DATA5 *)calloc(DIM+1, sizeof(DATA5));
	  d = d2;
	}
      
      // compute lik only if sel = 1
      nres = 0;
      maxL = 0;
      for(j=startT;j<=stopT;j++)
	{
	  if( Psel[j] == 1 ) 
	    {
	      for(ii=0;ii<NVAR;ii++)
		chi[ii] = 0;
	      
	      if(Pvar[0] == 1)
		chi[0] = (Teff[j] - wstar[0])/locsigma[0];
	      if(Pvar[1] == 1)
		chi[1] = (logg[j] - wstar[1])/locsigma[1];
	      if(Pvar[2] == 1)
		chi[2] = (z[j] - wstar[2])/locsigma[2];
	      if(Pvar[3] == 1)
		chi[3] = (Dni[j] - wstar[3])/locsigma[3];
	      if(Pvar[4] == 1)
		chi[4] = (nimax[j] - wstar[4])/locsigma[4];
	      if(Pvar[5] == 1)
		chi[5] = (M[j] - wstar[5])/locsigma[5];
	      if(Pvar[6] == 1)
		chi[6] = (R[j] - wstar[6])/locsigma[6];
	      
	      chi2 = 0;
	      for(n=0;n<NVAR;n++)
		chi2 += chi[n]*chi[n];
	      
	      // likelihood
	      L = mult * exp(-0.5*chi2);
	      if(L > maxL)
		maxL = L;
	      d[nres].L = L;
	      d[nres].M = M[j];
	      d[nres].R = R[j];
	      d[nres].logage = logage[j];
	      d[nres].pcage = pcage[j];
	      nres++;
	    }
	}
      if(nstar==0) 
	{
	  nres1 = nres;
	  maxL1 = maxL;
	}
      else
	{
	  nres2 = nres;
	  maxL2 = maxL;
	}
    }
  
   // independent estimates
  for(nstar=0;nstar<2;nstar++)
    {
      mass = radius = lt = ltnlog = rpcage = 0;
      count = 0;
      if(nstar==0)
	{
	  nres = nres1;
	  maxL = maxL1;
	  d = d1;
	}
      else
	{
	  nres = nres2;
	  maxL = maxL2;
	  d = d2;
	}
      
      // select only points with L >= 0.95 maxL
      for(j=0;j<nres;j++)
	{
	  if(d[j].L >= 0.95*maxL) 
	    {
	      mass += d[j].M;
	      radius += d[j].R;
	      lt += d[j].logage;
	      rpcage += d[j].pcage;
	      ltnlog += 1e-9*pow(10, d[j].logage);
	      count++;
	    }
	}
      mass /= (double)(count);
      radius /= (double)(count);
      lt /= (double)(count);
      ltnlog /= (double)(count);
      rpcage /= (double)(count);
     
      Pres[0+6*nstar] = mass;
      Pres[1+6*nstar] = radius;
      Pres[2+6*nstar] = lt;
      Pres[3+6*nstar] = ltnlog;
      Pres[4+6*nstar] = maxL;
      Pres[5+6*nstar] = rpcage;
    }

  // joint estimates
  qsort(d2, nres2, sizeof(DATA5), orderage);
  age2 = (double*)malloc(nres2*sizeof(double));
  age1 = (double*)malloc(nres1*sizeof(double));

  for(i=0;i<nres1;i++)
    age1[i] = 1e-9*pow(10, d1[i].logage);

  for(i=0;i<nres2;i++)
    age2[i] = 1e-9*pow(10, d2[i].logage);

  maxL = 0;
  for(j=0;j<nres1;j++)
    {
      findrange(age2, nres2, age1[j]-t_spread,age1[j]+t_spread, &lb, &ub);
      // the joint estimate is impossible
      if(lb == -1 || ub == -1) continue;
      if(lb == ub && fabs(age1[j] - age2[lb]) > t_spread) continue;
      for(i=lb;i<=ub;i++)
	{
	  count++;
	  L = d1[j].L * d2[i].L;
	  if(L > maxL)
	    maxL = L;
	}
    }

  for(j=12;j<22;j++)
    Pres[j] = 0;
  
  count = 0;
  for(j=0;j<nres1;j++)
    {
      findrange(age2, nres2, age1[j]-t_spread,age1[j]+t_spread, &lb, &ub);
      if(lb == -1 || ub == -1) continue;
      if(lb == ub && fabs(age1[j] - age2[lb]) > t_spread) continue;
      for(i=lb;i<=ub;i++)
	{
	  L = d1[j].L * d2[i].L;
	  if(L > 0.95*maxL)
	    {
	      
	      Pres[12] += d1[j].M;
	      Pres[13] += d1[j].R;
	      Pres[14] += d1[j].logage;
	      Pres[15] += age1[j];
	      
	      Pres[16] += d2[i].M;
	      Pres[17] += d2[i].R;
	      Pres[18] += d2[i].logage;
	      Pres[19] += age2[i];

	      Pres[21] += d1[j].pcage;

	      count++;
	    }
	}
    }
  
  Pres[20] = (double)count;

  for(j=12;j<20;j++)
    Pres[j] /= (double)(count);
  Pres[21] /= (double)(count);

  PROTECT( res = NEW_NUMERIC(22) );
  Rres = NUMERIC_POINTER(res);
  for(j=0;j<22;j++)
    Rres[j] = Pres[j];

  free(d1);
  free(d2);
  free(Psel);
  free(age1);
  free(age2);
  
  // exit
  UNPROTECT(9);
  return(res);
}
