/* wccsom.c: for mapping patters using the WCC criterion as a distance
   measure.
   This version: Jan 5, 2012
   Author: Ron Wehrens

   Adapted from:
   *
   *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-2002
   */

#include <R.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

double wcc_crosscorr(double *, double *, int, double *, int);
double wcc_autocorr(double *, int, double *, int);
double wcc_corr(double *, double *, int);


void WCC_onlineSOM(double *data, 
		   double *codes, 
		   double *nhbrdist,
		   double *alphas, double *pradius, 
		   int *ptrwdth, double *wghts,
		   double *dataAcors, double *Acors, 
		   double *changes,
		   int *pn, int *pp, 
		   int *pncodes, int *prlen)
{
  int n = *pn, p = *pp, ncodes = *pncodes, trwdth= *ptrwdth, rlen = *prlen; 
  int cd, i, j, k, l, nearest, fivers, twenties, iter;
  double dm, sim, dist, tmp, alpha, radius, decay;

  /* radius: fast exponential decay ending at zero; lateron we add
     1 so that at least one group of neighbours will be
     updated. Minimal value of 0.1 */
  if (*pradius > 1.1)
    decay = log(*pradius-1.0);
  else
    decay = 0.1;

  /* First calculate all autocors */
  for (i = 0; i < n; i++)
    dataAcors[i] = wcc_autocorr(data + i*p, p, wghts, trwdth);
  for (cd = 0; cd < ncodes; cd++)
    Acors[cd] = wcc_autocorr(codes + cd*p, p, wghts, trwdth);

  RANDIN;

  fivers = rlen*n / 20;
  twenties = rlen*n / 5;

  iter = 0;
  for (k = 0; k < rlen; k++) {
    alpha = alphas[1] + ((alphas[0] - alphas[1]) * (rlen - k) / rlen);
    radius = 1.0 + (*pradius-1.0) * exp(-decay * k / rlen);

    changes[k] = 0.0;

    for (l = 0; l < n; l++) {
      if (iter++ % twenties == 0)
	Rprintf("%d%%", 20 * iter / twenties);
      else if (iter % fivers == 0)
	Rprintf(".");
      R_FlushConsole();

      /* pick a random data point */
      i = (int)(n * UNIF);
      
      /* find the nearest code 'near' */
      nearest = -1;
      dm = -1.0;
      for (cd=0; cd<ncodes; cd++) {
	sim = wcc_crosscorr(data + i*p, 
			     codes + cd*p, 
			     p, wghts, trwdth)/(Acors[cd] * dataAcors[i]);
	if (sim > dm) {
	  nearest = cd;
	  dm = sim;
	}
      }

      if (nearest < 0)
	error("No nearest neighbour found...");

      /* update all codes within certain radius of 'nearest' */
      for (cd = 0; cd < ncodes; cd++) {
	if (nhbrdist[cd + ncodes*nearest] > radius) continue;
	tmp = alpha / (1.0 + nhbrdist[cd + ncodes*nearest]);

	for (j = 0; j < p; j++) {
	  dist = data[i*p + j] - codes[cd*p + j];
	  codes[cd*p + j] += tmp * dist;
	  
	  if (cd == nearest) changes[k] += dist * dist;
	}
	
	Acors[cd] = wcc_autocorr(codes + cd*p, p, wghts, trwdth);
	
	/*	if (cd == nearest) {
		tmp = wcc_crosscorr(data + i*p,
		codes + cd*p,
		p, wghts, trwdth)/(Acors[cd] * dataAcors[i]);
		changes[k] += tmp - dm; } */
      }
    }
  }
  
  for (k = 0; k < rlen; k++)
    changes[k] = sqrt(changes[k] / p)/n;
  /*    changes[k] /= n; */

  Rprintf("\n");
  R_FlushConsole();
  RANDOUT;
}


void WCCXYF_Tani(double *data, double *Ys, 
		 double *codes, double *codeYs,
		 double *nhbrdist,
		 double *alphas, double *pradius,	
		 double *xweight,
		 int *ptrwdth, double *wghts,
		 double *dataAcors, double *Acors,
		 double *changes,
		 double *xdists, double *ydists, /* working arrays */
		 int *pn, int *ppx, int *ppy, 
		 int *pncodes, int *prlen)
{
  int n = *pn, py = *ppy, px = *ppx, ncodes = *pncodes, rlen = *prlen,
    trwdth = *ptrwdth;
  int cd, i, j, k, l, nearest, fivers, twenties, iter;
  double dist, tmp, alpha, decay, radius, maxx, maxy;

  /* radius: exponential decay, after one-third of the iterations
     smaller than one */
  if (*pradius > 1.0)
    decay = 3.0 * log(*pradius);
  else
    decay = 1.0;

  /* First calculate all autocors */
  for (i = 0; i < n; i++)
    dataAcors[i] = wcc_autocorr(data + i*px, px, wghts, trwdth);
  for (cd = 0; cd < ncodes; cd++)
    Acors[cd] = wcc_autocorr(codes + cd*px, px, wghts, trwdth);

  RANDIN;

  fivers = rlen*n / 20;
  twenties = rlen*n / 5;
  iter = 0;

  for (k = 0; k < rlen; k++) {
    /* linear decrease in alpha, exponential in radius */
    alpha = alphas[1] + ((alphas[0] - alphas[1]) * (rlen - k) / rlen);
    radius = *pradius * exp(-decay * k / rlen);

    changes[k] = 0.0;

    for (l = 0; l < n; l++) {
      if (iter++ % twenties == 0)
	Rprintf("%d%%", 20 * iter / twenties);
      else if (iter % fivers == 0)
	Rprintf(".");
      R_FlushConsole();

      /* i is a counter over objects in data, cd is a counter over units
	 in the map, and j is a counter over variables */
      i = (int)(n * UNIF);
      
      /* calculate distances in x and y spaces. Both are scaled
	 between 0 and 1. The x-distance is a correlation measure
	 which is converted to a distance by taking 1-r. */
      maxx = maxy = 0;
      for (cd = 0; cd < ncodes; cd++) {
	xdists[cd] = 
	  1.0 - wcc_crosscorr(data + i*px, 
			      codes + cd*px, 
			      px, wghts, trwdth)/(Acors[cd] * dataAcors[i]);
	if (xdists[cd] > maxx) maxx = xdists[cd];
	
	/* Tanimoto distance */
	tmp = 0;
	for (j = 1; j < py; j++) {
	  if ((Ys[i*py + j] > .5 && 
	       codeYs[cd*py + j] <= .5) ||
	      (Ys[i*py + j] <= .5 && 
	       codeYs[cd*py + j] > .5))
	    tmp += 1.0;
	}
	ydists[cd] = tmp/(double)py;
	if (ydists[cd] > maxy) maxy = ydists[cd];
	/*	fprintf(stderr, "\nTanimoto distance: %.4lf",
		ydists[cd]); */
      }

      /* the following should never occur for x but may occur in freak
	 classification situations where all units end up in the same
	 class. */
      if (maxx < 1e-5) maxx = 1.0;
      if (maxy < 1e-5) maxy = 1.0;
      
      /* Find smallest distance. */
      dist = DOUBLE_XMAX; nearest = -1;
      for (cd = 0; cd < ncodes; cd++) {
	xdists[cd] /= maxx;
	ydists[cd] /= maxy;
	tmp = *xweight * xdists[cd] + (1.0 - *xweight) * ydists[cd];
	if (tmp < dist) {
	  dist = tmp;
	  nearest = cd;
	}
      }
      
      if (nearest < 0)
	error("No nearest neighbour found...");

     /* update all codes within certain radius of 'nearest' */
      for (cd = 0; cd < ncodes; cd++) {
	if(nhbrdist[cd + ncodes*nearest] > radius) continue;
	
	tmp = alpha / (1.0 + nhbrdist[cd + ncodes*nearest]);
	for(j = 0; j < px; j++) {
	  dist = data[i*px + j] - codes[cd*px + j];
	  codes[cd*px + j] += tmp * dist;
	  
	  if (cd == nearest) changes[k] += dist * dist;
	}
	
	/*	fprintf(stderr, "\n%d: Y\tcodeY", cd); */
	for(j = 0; j < py; j++) {
	  /*	  fprintf(stderr, "\n%.3lf\t%.3lf", Ys[i*py + j],
		  codeYs[cd*py + j]); */
	  dist = Ys[i*py + j] - codeYs[cd*py + j];
	  codeYs[cd*py + j] += tmp * dist;
	  
	  if (cd == nearest) changes[k+rlen] += dist * dist;
	}

	Acors[cd] = wcc_autocorr(codes + cd*px, px, wghts, trwdth);
      }
    }
  }
  
  for (k = 0; k < rlen; k++) {
    changes[k] = sqrt(changes[k]/px)/n;
    changes[k + rlen] = sqrt(changes[k + rlen]/py)/n;
  }

  Rprintf("\n");
  R_FlushConsole();
  RANDOUT;
}


void WCCXYF_Eucl(double *data, double *Ys, 
		 double *codes, double *codeYs,
		 double *nhbrdist,
		 double *alphas, double *pradius,	
		 double *xweight,
		 int *ptrwdth, double *wghts,
		 double *dataAcors, double *Acors,
		 double *changes,
		 double *xdists, double *ydists, /* working arrays */
		 int *pn, int *ppx, int *ppy, 
		 int *pncodes, int *prlen)
{
  int n = *pn, py = *ppy, px = *ppx, ncodes = *pncodes, rlen = *prlen,
    trwdth = *ptrwdth;
  int cd, i, j, k, l, nearest, fivers, twenties, iter;
  double dist, tmp, alpha, decay, radius, maxx, maxy;
  
  /* radius: exponential decay, after one-third of the iterations
     smaller than one */
  if (*pradius > 1.0)
    decay = 3.0 * log(*pradius);
  else
    decay = 1.0;
  
  /* First calculate all autocors */
  for (i = 0; i < n; i++)
    dataAcors[i] = wcc_autocorr(data + i*px, px, wghts, trwdth);
  for (cd = 0; cd < ncodes; cd++)
    Acors[cd] = wcc_autocorr(codes + cd*px, px, wghts, trwdth);
  
  RANDIN;
  
  fivers = rlen*n / 20;
  twenties = rlen*n / 5;
  iter = 0;
  
  for (k = 0; k < rlen; k++) {
    /* linear decrease in alpha, exponential in radius */
    alpha = alphas[1] + ((alphas[0] - alphas[1]) * (rlen - k) / rlen);
    radius = *pradius * exp(-decay * k / rlen);
    
    changes[k] = 0.0;
    
    for (l = 0; l < n; l++) {
      if (iter++ % twenties == 0)
	Rprintf("%d%%", 20 * iter / twenties);
      else if (iter % fivers == 0)
	Rprintf(".");
      R_FlushConsole();
      
      /* i is a counter over objects in data, cd is a counter over units
	 in the map, and j is a counter over variables */
      i = (int)(n * UNIF);
      
      /* calculate distances in x and y spaces. Both are scaled
	 between 0 and 1. The x-distance is a correlation measure
	 which is converted to a distance by taking 1-r. */
      maxx = maxy = -1.0;
      for (cd = 0; cd < ncodes; cd++) {
	xdists[cd] = 
	  1.0 - wcc_crosscorr(data + i*px, 
			      codes + cd*px, 
			      px, wghts, trwdth)/(Acors[cd] * dataAcors[i]);
	if (xdists[cd] > maxx) maxx = xdists[cd];
	
	/* Euclidean distance */
	ydists[cd] = 0.0;
	for (j = 0; j < py; j++) {
	  tmp = Ys[i*py + j] - codeYs[cd*py + j];
	  ydists[cd] += tmp * tmp;
	}
	ydists[cd] = sqrt(ydists[cd]);
	if (ydists[cd] > maxy) maxy = ydists[cd];
      }
      
      /* the following should never occur for x but may occur in freak
	 classification situations where all units end up in the same
	 class. */
      if (maxx < 1e-5) maxx = 1.0;
      if (maxy < 1e-5) maxy = 1.0;
      
      /* scaling of y distances in this case is necessary; we divide
	 by the largest distance. Then, add, with factor xweight, and
	 find smallest overall distance. */
      dist = DOUBLE_XMAX; nearest = -1;
      for (cd = 0; cd < ncodes; cd++) {
	xdists[cd] /= maxx;
	ydists[cd] /= maxy;
	tmp = *xweight * xdists[cd] + ((1.0 - *xweight) * ydists[cd]);
	if (tmp < dist) {
	  dist = tmp;
	  nearest = cd;
	}
      }
      
      if (nearest < 0)
	error("No nearest neighbour found...");

      /* update all codes within certain radius of 'nearest' */
      for (cd = 0; cd < ncodes; cd++) {
	if(nhbrdist[cd + ncodes*nearest] > radius) continue;
	
	tmp = alpha / (1.0 + nhbrdist[cd + ncodes*nearest]);
	for(j = 0; j < px; j++) {
	  dist = data[i*px + j] - codes[cd*px + j];
	  codes[cd*px + j] += tmp * dist;
	  
	  if (cd == nearest) changes[k] += dist * dist;
	}
	
	for(j = 0; j < py; j++) {
	  dist = Ys[i*py + j] - codeYs[cd*py + j];
	  codeYs[cd*py + j] += tmp * dist;
	  
	  if (cd == nearest) changes[k+rlen] += dist * dist;
	}

	Acors[cd] = wcc_autocorr(codes + cd*px, px, wghts, trwdth);
      }
    }
  }
  
  for (k = 0; k < rlen; k++) {
    changes[k] = sqrt(changes[k]/px)/n;
    changes[k + rlen] = sqrt(changes[k + rlen]/py)/n;
  }
  
  Rprintf("\n");
  R_FlushConsole();
  RANDOUT;
}



/* knn1 procedure using the wcc value as a distance */

void wccassign(double *data, double *dataAcors, 
	       double *codes, double *Acors, 
	       int *ptrwdth, double *wghts, 
	       int *classif, double *bestwccs,
	       int *pn, int *pp, int *pncodes)
{
  int n = *pn, p = *pp, ncodes = *pncodes, trwdth= *ptrwdth;
  int cd, i, nearest;
  double dm, wccval;

  for (i=0; i<n; i++) {
    nearest = -1;
    dm = 0.0;

    for (cd=0; cd<ncodes; cd++) {
      wccval = wcc_crosscorr(data + i*p, codes + cd*p, p, wghts, trwdth) / 
	(Acors[cd]*dataAcors[i]);

      if (wccval > dm) {
	nearest = cd;
	dm = wccval;
      }
    }

    classif[i] = nearest+1;
    bestwccs[i] = dm;
  }
}


double wcc_crosscorr(double *p1, double *p2, int np, 
		     double *wghts, int trwdth) {
  int i;
  double crosscov;

  crosscov = wcc_corr(p1, p2, np);
  
  for (i=1; i<trwdth; i++) {
    crosscov+=(wcc_corr(p1, p2+i, np-i)*wghts[i]);
    crosscov+=(wcc_corr(p1+i, p2, np-i)*wghts[i]);
  }
    
  return(crosscov);
}

double wcc_autocorr(double *p1, int np, double *wghts, int trwdth) {
  int i;
  double autocov;

  autocov = wcc_corr(p1, p1, np);

  for (i=1; i<trwdth; i++) 
    autocov += 2*(wcc_corr(p1, p1+i, np-i)*wghts[i]);

  return(sqrt(autocov));
}

double wcc_corr(double *p1, double *p2, int npoints)
{
  int i;
  double anum;

  anum = 0.0;
  for (i=0; i<npoints; i++)
    anum += p1[i]*p2[i];

  return(anum);
}

void wccdist(double *p1, double *p2, int *pnpoints, 
	     double *wghts, int *ptrwdth, double *WCC)
{
  int npoints = *pnpoints, trwdth=*ptrwdth;

  *WCC = wcc_crosscorr(p1, p2, npoints, wghts, trwdth);
}

void wacdist(double *p1, int *pnpoints, 
	     double *wghts, int *ptrwdth, double *ACC) 
{
  int npoints = *pnpoints, trwdth=*ptrwdth;

  *ACC = wcc_autocorr(p1, npoints, wghts, trwdth);
}


void wacdists(double *patterns, int *pnobj, int *pnpoints, 
	      double *wghts, int *ptrwdth, double *ACC) 
{
  int nobj = *pnobj, npoints = *pnpoints, trwdth=*ptrwdth;
  int i;

  for (i=0; i<nobj; i++)
    ACC[i] = wcc_autocorr(patterns + i*npoints, npoints, wghts, trwdth);
}

