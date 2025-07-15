#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>

void period();
void mset();
void eset();
void rset();
void spinset();
void sineset();
void mc();
void metro();
void single_clus();
void hybrid_update();
void init_histogram();
void update_histogram(double energy);
double analyze_histogram();
int detect_two_peaks();

extern void ranset(int, int []);
extern void rnd(int [], int, int []);

#define PI 3.141592653589793
#define iqq 8  // 8-state cube spins (vertices of a cube)
#define irbit 2281
#define mask 017777777777
#define update "hybrid"
#define spin_model "cube"  // Using cube spin model

int nx, ny, nla;  // Fixed at 3 layers

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
double fm_layer[3][8];  // thermodynamic quantities for each layer

// First-order transition analysis variables
#define HIST_SIZE 1000
int energy_hist[HIST_SIZE];        // Energy histogram
double e_min, e_max, de_hist;      // Energy histogram parameters
double latent_heat;                // Latent heat for first-order transition
int two_peak_detected;             // Flag for two-peak structure in histogram

// Hybrid Monte Carlo algorithm parameters
double metropolis_fraction;        // Fraction of updates that should be Metropolis
int cluster_interval;              // Interval between cluster updates

// Cube model parameters
int current_spin_model;             // Always 0 for cube model

int *ir;
int *irsd1;

int main() {
    int iri, i, itemp, layer;
    double fnla2, fnla4, temp;
    double fm2,fm4,fg2,fg4,fe1,fe2,cv;
    double ff,corr;
    int input_nx, input_ny;
    char dirname[256], filename[256];
    FILE *layer_files[3];

    // Read nx, ny, nmcs1, nmcs2, iri from input (fixed at 3 layers)
    scanf("%d %d %d %d %d", &nx, &ny, &nmcs1, &nmcs2, &iri);
    int nz = 3;  // Fixed at 3 layers
    nla = nx * ny * nz;  // Total sites = nx * ny * 3
    fnla2 = nla * nla;
    fnla4 = fnla2 * fnla2;

    printf("#%12d %12d %12d %12d %12d %12d\n", nx, ny, nz, nmcs1, nmcs2, iri);
    
    // Create timestamped directory structure
    time_t now = time(NULL);
    struct tm *local_time = localtime(&now);
    char timestamp[32];
    sprintf(timestamp, "%04d%02d%02d_%02d%02d%02d", 
            local_time->tm_year + 1900,
            local_time->tm_mon + 1,
            local_time->tm_mday,
            local_time->tm_hour,
            local_time->tm_min,
            local_time->tm_sec);
    
    // Create hierarchical directory structure: simulation_runs/timestamp/lattice_size/
    char base_dir[512], size_dir[512];
    sprintf(base_dir, "simulation_runs/%s_%s_%s", timestamp, spin_model, update);
    sprintf(size_dir, "%s/nx%d_ny%d", base_dir, nx, ny);
    sprintf(dirname, "%s", size_dir);
    
    mkdir("simulation_runs", 0755);
    mkdir(base_dir, 0755);
    mkdir(size_dir, 0755);

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

    // Open output files for each layer
    for(layer=0; layer<3; layer++) {
        sprintf(filename, "%s/layer_%d.txt", dirname, layer+1);
        layer_files[layer] = fopen(filename, "w");
        if (layer_files[layer] == NULL) {
            printf("Error: Cannot create file %s\n", filename);
            return 1;
        }
        // Write enhanced model metadata
        fprintf(layer_files[layer], "# Enhanced Quasi-3D Monte Carlo Simulation\n");
        fprintf(layer_files[layer], "# Timestamp: %s\n", timestamp);
        fprintf(layer_files[layer], "# Spin model: %s (%d states - cube vertices)\n", spin_model, iqq);
        fprintf(layer_files[layer], "# Algorithm: %s (hybrid Metropolis+Wolff)\n", update);
        fprintf(layer_files[layer], "# Theoretical Tc (8-state Potts): ~0.751\n");
        fprintf(layer_files[layer], "# Temperature range: 0.1-1.5 (optimized for critical region)\n");
        fprintf(layer_files[layer], "# Layer: %d/3\n", layer+1);
        fprintf(layer_files[layer], "# System size: %dx%dx3, Total sites: %d\n", nx, ny, nla);
        fprintf(layer_files[layer], "# Spin vectors: 8 cube vertices (±1,±1,±1)/√3\n");
        fprintf(layer_files[layer], "#\n");
        fprintf(layer_files[layer], "#%12s %12s %12s %12s %12s %12s %12s %12s\n",
                "Temperature", "M^2", "M^4", "G^2", "G^4", "Energy", "Cv", "Corr");
    }

    // Updated temperature range based on 8-state Potts theory
    // Theoretical Tc = 1/log(1+sqrt(8)) ≈ 0.751 for 8-state Potts model
    // Focus on critical region with finer resolution near Tc
    for(itemp=1; itemp<=201; itemp++)
    {
      // Temperature range: 0.1 to 1.5 (avoiding T=0 singularity)
      // Higher resolution near critical temperature
      if (itemp <= 50) {
        // Low temperature region: 0.1 to 0.5
        temp = 0.1 + 0.008*(itemp-1);
      } else if (itemp <= 150) {
        // Critical region: 0.5 to 1.0 (fine resolution around Tc ≈ 0.751)
        temp = 0.5 + 0.005*(itemp-50);
      } else {
        // High temperature region: 1.0 to 1.5
        temp = 1.0 + 0.01*(itemp-150);
      }
      
      beta=1/temp;
      
      // Initialize histogram for first-order transition analysis
      if (itemp == 1) {
        init_histogram();
      }
      
      spinset();

      mc();  // This now fills fm_layer[layer][i] for each layer

      // Calculate and output thermodynamic properties for each layer
      for(layer=0; layer<3; layer++) {
          double layer_size = nx * ny;  // Sites per layer
          double fnla2_layer = layer_size * layer_size;
          double fnla4_layer = fnla2_layer * fnla2_layer;
          
          fm2=fm_layer[layer][1]/fnla2_layer;
          fm4=fm_layer[layer][2]/fnla4_layer;
          fg2=fm_layer[layer][3]/layer_size;
          fg4=fm_layer[layer][4]/layer_size;
          fe1=fm_layer[layer][5]/layer_size;
          fe2=fm_layer[layer][6]/fnla2_layer;
          cv =beta*beta*(fe2-fe1*fe1)*layer_size;
          ff =fm_layer[layer][7]/fnla2_layer;
          
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
          
          fprintf(layer_files[layer], "%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                  temp,fm2,fm4,fg2,fg4,fe1,cv,corr);
      }
      
      // First-order transition analysis
      // Check for two-peak structure in energy distribution near critical temperature
      if (temp >= 0.6 && temp <= 0.9) {  // Near theoretical Tc ≈ 0.751
        two_peak_detected = detect_two_peaks();
        if (two_peak_detected) {
          latent_heat = analyze_histogram();
          printf("# First-order transition detected at T=%.4f, Latent heat=%.6f\n", 
                 temp, latent_heat);
        }
      }
    }
    
    // Close layer files
    for(layer=0; layer<3; layer++) {
        fclose(layer_files[layer]);
    }

    printf("Simulation completed. Layer data saved to %s/\n", dirname);
    
    // Free allocated memory
    free(isp); free(nn); free(n2); free(n4);
    free(cosx); free(sinx); free(cosy); free(siny);
    free(ir); free(irsd1);
    return 0;
}

void period()
/*
       periodic boundary conditions for quasi-3d system
       system size = nx*ny*nz (3 layers)
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
      nn[la+4*nla] = layer_pos + ((iz+1)%3)*layer_size;     // layer above
      nn[la+5*nla] = layer_pos + ((iz-1+3)%3)*layer_size;   // layer below
      
      // For n2 and n4 arrays (used for longer-range correlations)
      n2[la] = layer_pos + ((iz+1)%3)*layer_size;            // next layer in 3-layer system
      n4[la] = layer_pos + ((iz+2)%3)*layer_size;            // layer after next in 3-layer system
      n2[la+nla] = ((ix+nx/2)%nx) + iy*nx + iz*layer_size;   // half-way across x
      n4[la+nla] = ((ix+nx/4)%nx) + iy*nx + iz*layer_size;   // quarter-way across x
    }
}

void mset()
/*
        set magnetization vectors for 8-state cube spins
        8 vertices of a cube: (±1,±1,±1) normalized
*/
{
  int iq;
  double invsqrt3 = 1.0/sqrt(3.0);
  
  // Set spin model to cube
  current_spin_model = 0;
  
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
        rule for energy - 8-state cube model
*/
{
    int iq1,iq2;

    for (iq1=0; iq1 <= iqq-1; iq1++){
      for (iq2=0; iq2 <= iqq-1; iq2++){
        rule[iq1][iq2] = - (mx[iq1]*mx[iq2] + my[iq1]*my[iq2] + mz[iq1]*mz[iq2]);
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

  // Initialize hybrid algorithm parameters
  metropolis_fraction = 0.7;  // 70% Metropolis, 30% cluster updates
  cluster_interval = (int)(1.0 / (1.0 - metropolis_fraction));

  for (mcs=1; mcs <= nmcs1; mcs++){

      if(strcmp(update,"me")==0) { metro(); }
      if(strcmp(update,"wolff")==0) { single_clus(); }
      if(strcmp(update,"hybrid")==0) { hybrid_update(); }
  }

/*   measurement */

  for(i=1; i<=7; i++){
    fm[i] = 0;
    for(int layer=0; layer<3; layer++) {
      fm_layer[layer][i] = 0;
    }
  }

  for (mcs=1; mcs <= nmcs2; mcs++){

      if(strcmp(update,"me")==0) { metro(); }
      if(strcmp(update,"wolff")==0) { single_clus(); }
      if(strcmp(update,"hybrid")==0) { hybrid_update(); }

/*  measurement of order parameter, energy - layer by layer */

      // Initialize sums for each layer
      double fmxsum_layer[3], fmysum_layer[3], fmzsum_layer[3];
      double fm2xsum_layer[3], fm2ysum_layer[3], fm2zsum_layer[3];
      double fm4xsum_layer[3], fm4ysum_layer[3], fm4zsum_layer[3];
      
      for(int layer=0; layer<3; layer++) {
        fmxsum_layer[layer] = 0; fmysum_layer[layer] = 0; fmzsum_layer[layer] = 0;
        fm2xsum_layer[layer] = 0; fm2ysum_layer[layer] = 0; fm2zsum_layer[layer] = 0;
        fm4xsum_layer[layer] = 0; fm4ysum_layer[layer] = 0; fm4zsum_layer[layer] = 0;
      }
      
      // Overall system sums for backward compatibility
      fmxsum=0; fmysum=0; fmzsum=0;
      fm2xsum=0; fm2ysum=0; fm2zsum=0;
      fm4xsum=0; fm4ysum=0; fm4zsum=0;

      for (la=0; la <= nla-1; la++){
        int layer_size = nx * ny;
        int current_layer = la / layer_size;  // Determine which layer this site belongs to
        
        isp1=isp[la];
        
        // Add to layer-specific sums
        fmxsum_layer[current_layer] += mx[isp1];
        fmysum_layer[current_layer] += my[isp1];
        fmzsum_layer[current_layer] += mz[isp1];
        
        // Add to overall sums
        fmxsum += mx[isp1];
        fmysum += my[isp1];
        fmzsum += mz[isp1];
        
        // Calculate correlations within layer and between layers
        isp2 = isp[n2[la]];
        fm2xsum_layer[current_layer] += mx[isp1]*mx[isp2];
        fm2ysum_layer[current_layer] += my[isp1]*my[isp2];
        fm2zsum_layer[current_layer] += mz[isp1]*mz[isp2];
        fm2xsum += mx[isp1]*mx[isp2];
        fm2ysum += my[isp1]*my[isp2];
        fm2zsum += mz[isp1]*mz[isp2];
        
        isp2 = isp[n2[la+nla]];
        fm2xsum_layer[current_layer] += mx[isp1]*mx[isp2];
        fm2ysum_layer[current_layer] += my[isp1]*my[isp2];
        fm2zsum_layer[current_layer] += mz[isp1]*mz[isp2];
        fm2xsum += mx[isp1]*mx[isp2];
        fm2ysum += my[isp1]*my[isp2];
        fm2zsum += mz[isp1]*mz[isp2];
        
        isp2 = isp[n4[la]];
        fm4xsum_layer[current_layer] += mx[isp1]*mx[isp2];
        fm4ysum_layer[current_layer] += my[isp1]*my[isp2];
        fm4zsum_layer[current_layer] += mz[isp1]*mz[isp2];
        fm4xsum += mx[isp1]*mx[isp2];
        fm4ysum += my[isp1]*my[isp2];
        fm4zsum += mz[isp1]*mz[isp2];
        
        isp2 = isp[n4[la+nla]];
        fm4xsum_layer[current_layer] += mx[isp1]*mx[isp2];
        fm4ysum_layer[current_layer] += my[isp1]*my[isp2];
        fm4zsum_layer[current_layer] += mz[isp1]*mz[isp2];
        fm4xsum += mx[isp1]*mx[isp2];
        fm4ysum += my[isp1]*my[isp2];
        fm4zsum += mz[isp1]*mz[isp2];
      }
      
      // Calculate order parameters for each layer
      for(int layer=0; layer<3; layer++) {
        double f2order_layer = (fmxsum_layer[layer]*fmxsum_layer[layer] + 
                               fmysum_layer[layer]*fmysum_layer[layer] + 
                               fmzsum_layer[layer]*fmzsum_layer[layer]);
        double g2order_layer = (fm2xsum_layer[layer] + fm2ysum_layer[layer] + fm2zsum_layer[layer]);
        double g4order_layer = (fm4xsum_layer[layer] + fm4ysum_layer[layer] + fm4zsum_layer[layer]);
        
        fm_layer[layer][1] += f2order_layer;
        fm_layer[layer][2] += f2order_layer*f2order_layer;
        fm_layer[layer][3] += g2order_layer/2;
        fm_layer[layer][4] += g4order_layer/2;
      }
      
      // Overall system calculations for backward compatibility
      f2order=(fmxsum*fmxsum+fmysum*fmysum+fmzsum*fmzsum);
      g2order=(fm2xsum+fm2ysum+fm2zsum);
      g4order=(fm4xsum+fm4ysum+fm4zsum);

      fm[1] += f2order;
      fm[2] += f2order*f2order;
      fm[3] += g2order/2;
      fm[4] += g4order/2;

      // Energy calculation for each layer
      double fenergy_layer[3] = {0, 0, 0};
      fenergy=0;
      
      for (la=0; la <= nla-1; la++){
        int layer_size = nx * ny;
        int current_layer = la / layer_size;
        
        isp1=isp[la];
        double site_energy = rule[isp1][isp[nn[la]]]
                           + rule[isp1][isp[nn[la+2*nla]]]
                           + rule[isp1][isp[nn[la+4*nla]]]
                           + rule[isp1][isp[nn[la+5*nla]]];
        
        fenergy_layer[current_layer] += site_energy;
        fenergy += site_energy;
      }
      
      // Store energy data for each layer
      for(int layer=0; layer<3; layer++) {
        fm_layer[layer][5] += fenergy_layer[layer];
        fm_layer[layer][6] += fenergy_layer[layer]*fenergy_layer[layer];
      }
      
      fm[5] += fenergy;
      fm[6] += fenergy*fenergy;

      // Update energy histogram for first-order transition analysis
      update_histogram(fenergy);

      // Correlation calculations for each layer
      double clxc1_layer[3], clyc1_layer[3], clzc1_layer[3];
      double clxs1_layer[3], clys1_layer[3], clzs1_layer[3];
      double clxc2_layer[3], clyc2_layer[3], clzc2_layer[3];
      double clxs2_layer[3], clys2_layer[3], clzs2_layer[3];
      
      for(int layer=0; layer<3; layer++) {
        clxc1_layer[layer] = 0; clyc1_layer[layer] = 0; clzc1_layer[layer] = 0;
        clxs1_layer[layer] = 0; clys1_layer[layer] = 0; clzs1_layer[layer] = 0;
        clxc2_layer[layer] = 0; clyc2_layer[layer] = 0; clzc2_layer[layer] = 0;
        clxs2_layer[layer] = 0; clys2_layer[layer] = 0; clzs2_layer[layer] = 0;
      }
      
      clxc1 = 0; clyc1 = 0; clzc1 = 0;
      clxs1 = 0; clys1 = 0; clzs1 = 0;
      clxc2 = 0; clyc2 = 0; clzc2 = 0;
      clxs2 = 0; clys2 = 0; clzs2 = 0;

      for (la=0; la <= nla-1; la++){
        int layer_size = nx * ny;
        int current_layer = la / layer_size;
        int layer_pos = la % layer_size;
        
        isp1 = isp[la];
        
        // Layer-specific calculations
        clxc1_layer[current_layer] += mx[isp1]*cosx[layer_pos%nx];
        clyc1_layer[current_layer] += my[isp1]*cosx[layer_pos%nx];
        clzc1_layer[current_layer] += mz[isp1]*cosx[layer_pos%nx];
        clxs1_layer[current_layer] += mx[isp1]*sinx[layer_pos%nx];
        clys1_layer[current_layer] += my[isp1]*sinx[layer_pos%nx];
        clzs1_layer[current_layer] += mz[isp1]*sinx[layer_pos%nx];
        clxc2_layer[current_layer] += mx[isp1]*cosx[layer_pos/nx];
        clyc2_layer[current_layer] += my[isp1]*cosx[layer_pos/nx];
        clzc2_layer[current_layer] += mz[isp1]*cosx[layer_pos/nx];
        clxs2_layer[current_layer] += mx[isp1]*sinx[layer_pos/nx];
        clys2_layer[current_layer] += my[isp1]*sinx[layer_pos/nx];
        clzs2_layer[current_layer] += mz[isp1]*sinx[layer_pos/nx];
        
        // Overall system calculations
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
      
      // Calculate correlation function for each layer
      for(int layer=0; layer<3; layer++) {
        double cl_layer = clxc1_layer[layer]*clxc1_layer[layer] + clyc1_layer[layer]*clyc1_layer[layer] + clzc1_layer[layer]*clzc1_layer[layer]
                        + clxs1_layer[layer]*clxs1_layer[layer] + clys1_layer[layer]*clys1_layer[layer] + clzs1_layer[layer]*clzs1_layer[layer]
                        + clxc2_layer[layer]*clxc2_layer[layer] + clyc2_layer[layer]*clyc2_layer[layer] + clzc2_layer[layer]*clzc2_layer[layer]
                        + clxs2_layer[layer]*clxs2_layer[layer] + clys2_layer[layer]*clys2_layer[layer] + clzs2_layer[layer]*clzs2_layer[layer];
        cl_layer /= 4;
        fm_layer[layer][7] += cl_layer;
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
    for(int layer=0; layer<3; layer++) {
      fm_layer[layer][i] /= nmcs2;
    }
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

// End of Wolff cluster update

void hybrid_update()
/*
   Hybrid Monte Carlo update combining Metropolis and Wolff cluster algorithms
   Based on research by Hasenbusch (2020) on icosahedral models
*/
{
  static int update_counter = 0;
  update_counter++;
  
  // Decide whether to use Metropolis or cluster update
  // Use a deterministic pattern to ensure proper mixing
  if (update_counter % cluster_interval == 0) {
    // Cluster update (less frequent, more global changes)
    single_clus();
  } else {
    // Metropolis update (more frequent, local changes)
    // Perform multiple local updates to match cluster update efficiency
    int local_updates = nla / 10;  // Number of local updates per step
    for (int i = 0; i < local_updates; i++) {
      metro();
    }
  }
}

// First-order transition analysis functions

void init_histogram()
/*
   Initialize energy histogram for first-order transition detection
*/
{
  int i;
  
  // Initialize histogram
  for (i = 0; i < HIST_SIZE; i++) {
    energy_hist[i] = 0;
  }
  
  // Set energy range for histogram (rough estimates)
  e_min = -3.0 * nla;  // Minimum possible energy
  e_max = 0.0;         // Maximum possible energy (all spins random)
  de_hist = (e_max - e_min) / HIST_SIZE;
  
  two_peak_detected = 0;
  latent_heat = 0.0;
}

void update_histogram(double energy)
/*
   Update energy histogram with current energy value
*/
{
  int bin;
  
  if (de_hist > 0) {
    bin = (int)((energy - e_min) / de_hist);
    if (bin >= 0 && bin < HIST_SIZE) {
      energy_hist[bin]++;
    }
  }
}

int detect_two_peaks()
/*
   Detect two-peak structure in energy histogram (signature of first-order transition)
   Returns 1 if two peaks detected, 0 otherwise
*/
{
  int i, peak_count = 0;
  int local_maxima[10];  // Store positions of local maxima
  int threshold = nmcs2 / 100;  // Minimum height for a peak
  
  // Find local maxima
  for (i = 1; i < HIST_SIZE - 1; i++) {
    if (energy_hist[i] > energy_hist[i-1] && 
        energy_hist[i] > energy_hist[i+1] && 
        energy_hist[i] > threshold) {
      if (peak_count < 10) {
        local_maxima[peak_count] = i;
        peak_count++;
      }
    }
  }
  
  // Check if we have at least two significant peaks
  if (peak_count >= 2) {
    // Additional check: peaks should be separated by a minimum distance
    int separation = local_maxima[1] - local_maxima[0];
    if (separation > HIST_SIZE / 10) {  // Peaks should be well separated
      return 1;
    }
  }
  
  return 0;
}

double analyze_histogram()
/*
   Analyze energy histogram to extract latent heat
   Returns latent heat value
*/
{
  int i, peak1 = -1, peak2 = -1;
  int max1 = 0, max2 = 0;
  
  // Find two highest peaks
  for (i = 0; i < HIST_SIZE; i++) {
    if (energy_hist[i] > max1) {
      max2 = max1;
      peak2 = peak1;
      max1 = energy_hist[i];
      peak1 = i;
    } else if (energy_hist[i] > max2) {
      max2 = energy_hist[i];
      peak2 = i;
    }
  }
  
  // Calculate latent heat as energy difference between peaks
  if (peak1 >= 0 && peak2 >= 0) {
    double e1 = e_min + peak1 * de_hist;
    double e2 = e_min + peak2 * de_hist;
    return fabs(e2 - e1) / nla;  // Per site latent heat
  }
  
  return 0.0;
}
