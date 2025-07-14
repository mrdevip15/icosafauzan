#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void period();
void mset();
void eset();
void rset();
void spinset();
void sineset();
void mc();
void metro();
void single_clus();

extern void ranset(int, int []);
extern void rnd(int [], int, int []);

#define PI 3.141592653589793
#define iqq 8
#define irbit 2281
#define mask 017777777777
#define update "wolff"

int nx, ny, nla;  // Fixed at 4 layers

int *isp;
int *nn;
int *n2;
int *n4;
double mx[iqq], my[iqq], mz[iqq];
double rule[iqq][iqq];
int rr[iqq][15];
int nmcs1, nmcs2;
double *cosx, *sinx, *cosy, *siny;
double beta;
double fm[8];

int *ir;
int *irsd1;

int main() {
    int iri, i, itemp;
    double fnla2, fnla4, temp;
    double fm2,fm4,fg2,fg4,fe1,fe2,cv;
    double ff,corr;
    int input_nx, input_ny;

    // Read nx, ny, nmcs1, nmcs2, iri from input (fixed at 4 layers)
    scanf("%d %d %d %d %d", &nx, &ny, &nmcs1, &nmcs2, &iri);
    int nz = 4;  // Fixed at 4 layers
    nla = nx * ny * nz;  // Total sites = nx * ny * 4
    fnla2 = nla * nla;
    fnla4 = fnla2 * fnla2;

    printf("#%12d %12d %12d %12d %12d %12d\n", nx, ny, nz, nmcs1, nmcs2, iri);

    // Dynamically allocate arrays
    isp = (int*)malloc(nla * sizeof(int));
    nn = (int*)malloc(6 * nla * sizeof(int));  // 6 neighbors: 4 in-plane + 2 inter-layer
    n2 = (int*)malloc(2 * nla * sizeof(int));
    n4 = (int*)malloc(2 * nla * sizeof(int));
    cosx = (double*)malloc(nx * sizeof(double));
    sinx = (double*)malloc(nx * sizeof(double));
    cosy = (double*)malloc(ny * sizeof(double));
    siny = (double*)malloc(ny * sizeof(double));
    ir = (int*)malloc(4 * nla * sizeof(int));
    irsd1 = (int*)malloc(irbit * sizeof(int));

    ranset(iri, irsd1);
    for(i = 0; i<=irbit-1; i++){
      irsd1[i] &= mask;
    }

    period();
    mset();
    eset();
    rset();
    sineset();

    for(itemp=1; itemp<=201; itemp++)
    {
      temp=0.0+0.01*(itemp-1);
      beta=1/temp;
      spinset();

      mc();

      fm2=fm[1]/fnla2;
      fm4=fm[2]/fnla4;
      fg2=fm[3]/nla;
      fg4=fm[4]/nla;
      fe1=fm[5]/nla;
      fe2=fm[6]/fnla2;
      cv =beta*beta*(fe2-fe1*fe1)*nla;
      ff =fm[7]/fnla2;
      
      // Check for NaN or infinite values
      if (isnan(fm2) || isinf(fm2)) fm2 = 0.0;
      if (isnan(fm4) || isinf(fm4)) fm4 = 0.0;
      if (isnan(fg2) || isinf(fg2)) fg2 = 0.0;
      if (isnan(fg4) || isinf(fg4)) fg4 = 0.0;
      if (isnan(fe1) || isinf(fe1)) fe1 = 0.0;
      if (isnan(cv) || isinf(cv)) cv = 0.0;
      if (isnan(ff) || isinf(ff)) ff = 0.0;
      
      // Safe correlation length calculation
      double corr_arg = fm2/ff - 1.0;
      double sin_term = 2.0*sin(PI/nx);
      
      if (ff > 1e-10 && corr_arg > 0.0 && sin_term > 1e-10) {
        corr = sqrt(corr_arg) / sin_term;
      } else {
        corr = 0.0;  // Set to zero if calculation would be invalid
      }
      printf("%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
              temp,fm2,fm4,fg2,fg4,fe1,cv,corr);
    }

    // Free allocated memory
    free(isp); free(nn); free(n2); free(n4);
    free(cosx); free(sinx); free(cosy); free(siny);
    free(ir); free(irsd1);
    return 0;
}

void period()
/*
       periodic boundary conditions for quasi-3d system
       system size = nx*ny*nz (4 layers)
       Each layer is 2D with periodic boundary conditions
       Inter-layer connections with periodic boundary in z-direction
*/
{
    int la, ix, iy, iz, layer_size;
    
    layer_size = nx * ny;  // Size of each 2D layer
    
    for (la=0; la <= nla-1; la++){
      // Calculate 3D coordinates
      iz = la / layer_size;           // z-coordinate (layer number)
      int layer_pos = la % layer_size; // position within the layer
      iy = layer_pos / nx;            // y-coordinate within layer
      ix = layer_pos % nx;            // x-coordinate within layer
      
      // In-plane neighbors (4 neighbors within the same layer)
      nn[la]       = ((ix+1)%nx) + iy*nx + iz*layer_size;    // right neighbor
      nn[la+nla]   = ((ix-1+nx)%nx) + iy*nx + iz*layer_size; // left neighbor
      nn[la+2*nla] = ix + ((iy+1)%ny)*nx + iz*layer_size;    // up neighbor
      nn[la+3*nla] = ix + ((iy-1+ny)%ny)*nx + iz*layer_size; // down neighbor
      
      // Inter-layer neighbors (2 neighbors in adjacent layers)
      nn[la+4*nla] = layer_pos + ((iz+1)%4)*layer_size;     // layer above
      nn[la+5*nla] = layer_pos + ((iz-1+4)%4)*layer_size;   // layer below
      
      // For n2 and n4 arrays (used for longer-range correlations)
      n2[la] = layer_pos + ((iz+2)%4)*layer_size;            // half-way across 4 layers (2 layers away)
      n4[la] = layer_pos + ((iz+1)%4)*layer_size;            // quarter-way across 4 layers (1 layer away)
      n2[la+nla] = ((ix+nx/2)%nx) + iy*nx + iz*layer_size;   // half-way across x
      n4[la+nla] = ((ix+nx/4)%nx) + iy*nx + iz*layer_size;   // quarter-way across x
    }
}

void mset()
/*
        set magnetization
        cube
*/
{
  int iq;
  double invsqrt3 = 1.0/sqrt(3.0);

  /* 8 vertices of a cube: (±1,±1,±1) normalized */
  mx[0] = invsqrt3;  my[0] = invsqrt3;  mz[0] = invsqrt3;   /* (1,1,1) */
  mx[1] = -invsqrt3; my[1] = invsqrt3;  mz[1] = invsqrt3;   /* (-1,1,1) */
  mx[2] = -invsqrt3; my[2] = -invsqrt3; mz[2] = invsqrt3;   /* (-1,-1,1) */
  mx[3] = invsqrt3;  my[3] = -invsqrt3; mz[3] = invsqrt3;   /* (1,-1,1) */
  mx[4] = invsqrt3;  my[4] = invsqrt3;  mz[4] = -invsqrt3;  /* (1,1,-1) */
  mx[5] = -invsqrt3; my[5] = invsqrt3;  mz[5] = -invsqrt3;  /* (-1,1,-1) */
  mx[6] = -invsqrt3; my[6] = -invsqrt3; mz[6] = -invsqrt3;  /* (-1,-1,-1) */
  mx[7] = invsqrt3;  my[7] = -invsqrt3; mz[7] = -invsqrt3;  /* (1,-1,-1) */
}

void eset()
/*
        rule for energy
*/
{
    int iq1,iq2;

    for (iq1=0; iq1 <= iqq-1; iq1++){
      for (iq2=0; iq2 <= iqq-1; iq2++){
        rule[iq1][iq2] = - (mx[iq1]*mx[iq2] + my[iq1]*my[iq2]
                          + mz[iq1]*mz[iq2]);
      }
    }
}

void rset()
/*
        reflection of spins for cube
*/
{
    int iq1, iq2, rchoice;
    int i1[15],i2[15];
    double a, b, c, w, wx, wy, wz;
    double x1, y1, z1, x2, y2, z2;

    /* Define 15 reflection operations for cube */
    /* Face reflections (6 operations) */
    i1[0]=0, i2[0]=1;   /* x-reflection */
    i1[1]=2, i2[1]=3;   /* x-reflection */
    i1[2]=4, i2[2]=5;   /* x-reflection */
    i1[3]=6, i2[3]=7;   /* x-reflection */
    
    i1[4]=0, i2[4]=3;   /* y-reflection */
    i1[5]=1, i2[5]=2;   /* y-reflection */
    i1[6]=4, i2[6]=7;   /* y-reflection */
    i1[7]=5, i2[7]=6;   /* y-reflection */
    
    i1[8]=0, i2[8]=4;   /* z-reflection */
    i1[9]=1, i2[9]=5;   /* z-reflection */
    i1[10]=2, i2[10]=6; /* z-reflection */
    i1[11]=3, i2[11]=7; /* z-reflection */
    
    /* Diagonal reflections (3 operations) */
    i1[12]=0, i2[12]=6; /* diagonal reflection */
    i1[13]=1, i2[13]=7; /* diagonal reflection */
    i1[14]=2, i2[14]=4; /* diagonal reflection */

    for (rchoice=0; rchoice <= 14; rchoice++){
      x1=mx[i1[rchoice]], y1=my[i1[rchoice]], z1=mz[i1[rchoice]];
      x2=mx[i2[rchoice]], y2=my[i2[rchoice]], z2=mz[i2[rchoice]];

      a=y1*z2-y2*z1;
      b=z1*x2-z2*x1;
      c=x1*y2-x2*y1;

      for (iq1=0; iq1 <= iqq-1; iq1++){
        w = a*mx[iq1]+b*my[iq1]+c*mz[iq1];
        if(w>-0.001 && w<0.001){ 
          rr[iq1][rchoice] = iq1;
//            printf("%d %d %d\n",rchoice,iq1,iq1);
        } else {
          for (iq2=0; iq2 <= iqq-1; iq2++){
            w = a*(mx[iq1]+mx[iq2])+b*(my[iq1]+my[iq2])+c*(mz[iq1]+mz[iq2]);
            wx=(mx[iq2]-mx[iq1])*b-(my[iq2]-my[iq1])*a;
            wy=(my[iq2]-my[iq1])*c-(mz[iq2]-mz[iq1])*b;
            wz=(mz[iq2]-mz[iq1])*a-(mx[iq2]-mx[iq1])*c;
            if(wx>-0.001 && wx<0.001 && wy>-0.001 && wy<0.001 
              && wz>-0.001 && wz<0.001 && w>-0.001 && w<0.001) { 
              rr[iq1][rchoice] = iq2;
//              printf("%d %d %d\n",rchoice,iq1,iq2);
            }
          }
        }
      }
    }

}

void spinset()
/*
        set initial spins
*/
{
    int la;

    for (la=0; la <= nla-1; la++){
      isp[la]=0;
    }
}

void sineset()
/* set sine function */
{
    int la;

    for (la=0; la <= nx-1; la++){
        cosx[la] = cos(2*PI*la/nx);
        sinx[la] = sin(2*PI*la/nx);
    }
}

void mc()
/*
        monte carlo update
*/
{
  int mcs, i, iq;
  double fmxsum, fmysum, fmzsum, fenergy;
  double fm2xsum, fm2ysum, fm2zsum, fm4xsum, fm4ysum, fm4zsum;
  int la, la1, isp1, isp2;
  double ener, f2order, g2order, g4order;
  double clxc1, clyc1, clzc1, clxs1, clys1, clzs1;
  double clxc2, clyc2, clzc2, clxs2, clys2, clzs2;
  double cl;

/*   initialization  */

  for (mcs=1; mcs <= nmcs1; mcs++){

      if(update=="me") { metro(); }
      if(update=="wolff") { single_clus(); }
  }

/*   measurement */

  for(i=1; i<=7; i++){
    fm[i] = 0;
  }

  for (mcs=1; mcs <= nmcs2; mcs++){

      if(update=="me") { metro(); }
      if(update=="wolff") { single_clus(); }

/*  measurement of order parameter, energy */

      fmxsum=0;
      fmysum=0;
      fmzsum=0;
      fm2xsum=0;
      fm2ysum=0;
      fm2zsum=0;
      fm4xsum=0;
      fm4ysum=0;
      fm4zsum=0;

      for (la=0; la <= nla-1; la++){
        isp1=isp[la];
        fmxsum += mx[isp1];
        fmysum += my[isp1];
        fmzsum += mz[isp1];
        isp2 = isp[n2[la]];
        fm2xsum += mx[isp1]*mx[isp2];
        fm2ysum += my[isp1]*my[isp2];
        fm2zsum += mz[isp1]*mz[isp2];
        isp2 = isp[n2[la+nla]];
        fm2xsum += mx[isp1]*mx[isp2];
        fm2ysum += my[isp1]*my[isp2];
        fm2zsum += mz[isp1]*mz[isp2];
        isp2 = isp[n4[la]];
        fm4xsum += mx[isp1]*mx[isp2];
        fm4ysum += my[isp1]*my[isp2];
        fm4zsum += mz[isp1]*mz[isp2];
        isp2 = isp[n4[la+nla]];
        fm4xsum += mx[isp1]*mx[isp2];
        fm4ysum += my[isp1]*my[isp2];
        fm4zsum += mz[isp1]*mz[isp2];
      }
      f2order=(fmxsum*fmxsum+fmysum*fmysum+fmzsum*fmzsum);
      g2order=(fm2xsum+fm2ysum+fm2zsum);
      g4order=(fm4xsum+fm4ysum+fm4zsum);

      fm[1] += f2order;
      fm[2] += f2order*f2order;
      fm[3] += g2order/2;
      fm[4] += g4order/2;

      fenergy=0;
      for (la=0; la <= nla-1; la++){
        isp1=isp[la];
        fenergy += rule[isp1][isp[nn[la]]]
                 + rule[isp1][isp[nn[la+2*nla]]]
                 + rule[isp1][isp[nn[la+4*nla]]]
                 + rule[isp1][isp[nn[la+5*nla]]];
      }
      fm[5] += fenergy;
      fm[6] += fenergy*fenergy;

      clxc1 = 0;
      clyc1 = 0;
      clzc1 = 0;
      clxs1 = 0;
      clys1 = 0;
      clzs1 = 0;
      clxc2 = 0;
      clyc2 = 0;
      clzc2 = 0;
      clxs2 = 0;
      clys2 = 0;
      clzs2 = 0;

      for (la=0; la <= nla-1; la++){
        isp1 = isp[la];
        clxc1 += mx[isp1]*cosx[la%nx];
        clyc1 += my[isp1]*cosx[la%nx];
        clzc1 += mz[isp1]*cosx[la%nx];
        clxs1 += mx[isp1]*sinx[la%nx];
        clys1 += my[isp1]*sinx[la%nx];
        clzs1 += mz[isp1]*sinx[la%nx];
        clxc2 += mx[isp1]*cosx[la/nx];
        clyc2 += my[isp1]*cosx[la/nx];
        clzc2 += mz[isp1]*cosx[la/nx];
        clxs2 += mx[isp1]*sinx[la/nx];
        clys2 += my[isp1]*sinx[la/nx];
        clzs2 += mz[isp1]*sinx[la/nx];
      }

      cl = clxc1*clxc1+clyc1*clyc1+clzc1*clzc1
          +clxs1*clxs1+clys1*clys1+clzs1*clzs1
          +clxc2*clxc2+clyc2*clyc2+clzc2*clzc2
          +clxs2*clxs2+clys2*clys2+clzs2*clzs2;
      cl /= 4;

      fm[7] += cl;
  }


  for(i=1; i<=7; i++){
    fm[i] /= nmcs2;
  }
}

void metro()
/*   Metropilis */
{
  int la, la1, isp1, isp2;
  double ener1, ener2;

      rnd(ir, 4*nla, irsd1);
  for(la1=0; la1 <= nla-1; la1++){
      la=ir[la1]%nla;
      isp1=isp[la];
      ener1 = rule[isp1][isp[nn[la]]]
              +rule[isp1][isp[nn[la+nla]]]
              +rule[isp1][isp[nn[la+2*nla]]]
              +rule[isp1][isp[nn[la+3*nla]]]
              +rule[isp1][isp[nn[la+4*nla]]]
              +rule[isp1][isp[nn[la+5*nla]]];
      isp2=ir[la1+nla]%iqq;
      ener2 = rule[isp2][isp[nn[la]]]
              +rule[isp2][isp[nn[la+nla]]]
              +rule[isp2][isp[nn[la+2*nla]]]
              +rule[isp2][isp[nn[la+3*nla]]]
              +rule[isp2][isp[nn[la+4*nla]]]
              +rule[isp2][isp[nn[la+5*nla]]];
      if(ener2-ener1<0){
        isp[la]=isp2;
      } else {
        if(exp(-beta*(ener2-ener1))*mask > ir[la1+2*nla]){
          isp[la]=isp2;
        }
      } 
  }
}


void single_clus()
/*   Wolff */
{
  int la, la1, i, in, ic, rchoice, isp1;
  double edif;
  int mark[nla], next[nla];
  double boltz;

  rnd(ir,4*nla,irsd1);
  rchoice = ir[nla]%15;

// clear marked sites of previous cluster

  for(la=0; la<=nla-1; la++)
  {
    mark[la]=0;
    next[la]=la;
  }

// start growing anew cluster around site la, picked at random

  la = ir[2*nla]%nla;
  ic = 0;
  in = ic;
  next[in] = la;
  mark[la] = 1;

  while(in<=ic)
  {
    la = next[in];                    // current "boundary" spin
    isp1 = isp[la];
    isp[la] = rr[isp1][rchoice];

    for(i=0; i<=5; i++)              // test all 6 neighbors (4 in-plane + 2 inter-layer)
    {
      la1 = nn[la+i*nla];
      if(mark[la1]!=1){                // already member of a cluster
        edif = rule[isp[la1]][isp[la]]-rule[isp[la1]][isp1];
        if(edif > 0){                  // spins anti-parallel
          boltz = exp(-beta*edif);
          if(boltz*mask <= ir[la+i*nla]){    //inactive,try next neighbor

// otherwise the bond is active and site la1 is added to the current culster

            ic ++;
            mark[la1] = 1;       // active, then new member of cluster
            next[ic] = la1;      // store location of new cluster member
          }
        }
      }                         // end of loop over neighbors
    }

// now all neighbors of spin at site la are checked

    in ++;
  }                      // cluster still contains spins whose
 }                        // neighbors have not been tested
