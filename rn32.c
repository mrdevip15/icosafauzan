#define  nbit 32
#define  p    2281
#define  q    1029
int r = p-q;

#include <stdlib.h>

void ranset(int init, int irseed[])
/*
       make 32bit random number seed
*/
{
      int i,j,pp,qq,jtp;

      int *iy = (int *)malloc(p * sizeof(int));

      pp = p - 2;
      qq = p - 2+ q;

      for(i=0; i<=p-1; i++){
        irseed[i] = 0;
      }

      for(i=0; i<=30; i++){
        iy[i] = (init >> i) & 1;
      }

      for(i=31; i<=p-1; i++){
        iy[i] = iy[i-31] ^ iy[i-13];
      }

      for(i=0; i<=p*nbit-1; i++){
        pp = (pp+1) % (p-1);
        qq = (qq+1) % (p-1);

        jtp = iy[pp] ^ iy[qq];
        iy[pp] = jtp;

      }

      for(i=0; i<=p-1; i++){
        for(j=0; j<=nbit-1; j++){
          pp = (pp+1) % (p-1);
          qq = (qq+1) % (p-1);
          jtp = iy[pp] ^ iy[qq];
          iy[pp] = jtp;
          irseed[i] = (irseed[i] << 1) | jtp;
        }
      }

      pp = p - 2;
      qq = p - 2+ q;

      for(i=1; i<=100000; i++){
          pp = (pp+1) % (p-1);
          qq = (qq+1) % (p-1);
          irseed[pp] = irseed[pp] ^ irseed[qq];
      }

      free(iy);
}

void rnd(int ir[], int n,int irseed[])
{
      int *iy = (int *)malloc(r * sizeof(int));
      int *iz = (int *)malloc(p * sizeof(int));

      int i,j,repeat,remain;

      repeat = n/r;
      remain = n - repeat*r;
      for(j=0; j<=repeat-1; j++){
        for(i=0; i<=r-1; i++){
          iy[i]=irseed[i] ^ irseed[i+q];
          ir[j*r+i]=iy[i];
        }

        for(i=0; i<=q-1; i++){
          irseed[i]=irseed[i+r];
        }
        for(i=0; i<=r-1; i++){
          irseed[i+q]=iy[i];
        }
      }

      if(remain != 0) {
        for(i=0; i<=remain-1; i++){
          iy[i]=irseed[i] ^ irseed[i+q];
          ir[repeat*r+i]=iy[i];
        }
        for(i=0; i<=p-remain-1; i++){
          iz[i]=irseed[i+remain];
        }
        for(i=0; i<=p-remain-1; i++){
          irseed[i]=iz[i];
        }
        for(i=0; i<=remain-1; i++){
          irseed[i+p-remain]=iy[i];
        }
      }

      free(iy);
      free(iz);
}

