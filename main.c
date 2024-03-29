#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

void gaussian_nbr( double guass_var[2]);
double get_acc(double omega_0, double x, double m, double eta, double r_v, double v);
void save_x_to_disk(double ** x_pos, int nbr_of_simulations, double timestep, int nbr_of_particles, double * x_mean, double * x_variance);
void save_v_to_disk(double ** v_particles, int nbr_of_simulations, double timestep, int nbr_of_particles, double * v_mean, double * v_variance);
void get_histogram_values(int j, int i, double tmp_v, double tmp_x);
void save_histogram_values_to_disk(int avg_loops);

double* v1; double* v2; double* v3; double* v4;
double* x1; double* x2; double* x3; double* x4;
int counter1, counter2, counter3, counter4;

double tau;
double v_th;
double m, temperature, k_b, c_0, eta;
double timestep;
  
int main() {


  srand(time(NULL));
  // GSL INITIALIZATION
  const gsl_rng_type *T;
  gsl_rng *q;
  // Initializations
  gsl_rng_env_setup();
  T = gsl_rng_default;
  q = gsl_rng_alloc(T);
  gsl_rng_set(q,time(NULL));

  // Declaring variables
  int nbr_of_simulations;
  double gauss_rv[2];
  double omega_0;
  double radius;
  double density;
  int nbr_of_particles;
  int avg_loops;
  double start_v;
  double start_x;
  double tmp_v;
  double tmp_a;
  double tmp_x;
  double r_v;
  double uni_var1, uni_var2, uni_var3, uni_var4;
  double sigma; 
  
  // Initializing variable
  nbr_of_simulations = 25000;
  timestep = 5*1.0E-6; // 1e5 and 500 simulations is good for v
  nbr_of_particles = 5;
  avg_loops = 500;
  radius = (2.79*1E-6)/2.0;
  density = 2.65 * 1E3;
  m = 4.0*M_PI*radius*radius*radius*density/3.0; // assuming particle is sphere-ish
  temperature = 297;
  tau = 48.5*1E-6;  // A
  //  tau = 147.3*1E-6; // B
  eta = 1.0/tau;
  c_0 = exp(-eta*timestep);
  k_b = 1.380*1E-23;
  omega_0 = 3000.0*2.0*M_PI;
  v_th = sqrt(k_b*temperature / m);
  start_v = 2*1E-3;
  start_x = 1*1E-7;
  counter1=0;
  counter2=0;
  counter3=0;
  counter4=0;
  sigma = sqrt(2* eta * k_b * temperature / m);

  FILE * scale;
  scale = fopen("scale.dat", "w");
  fprintf(scale, "%f", timestep*1E3);
  fclose(scale);

  // Allocation space to save variables heap and stack
  double ** x_pos = malloc(sizeof(double)*(nbr_of_simulations));
  double * x_pointer = malloc(sizeof(double)*nbr_of_simulations*nbr_of_particles);
  double ** v_particles = malloc(sizeof(double)*(nbr_of_simulations));
  double * v_pointer = malloc(sizeof(double)*nbr_of_simulations*nbr_of_particles);
  double * v_mean = calloc(sizeof(double),nbr_of_simulations);
  double * v_variance = calloc(sizeof(double),nbr_of_simulations);
  double * x_mean = calloc(sizeof(double),nbr_of_simulations);
  double * x_variance = calloc(sizeof(double),nbr_of_simulations);

  double x[nbr_of_particles];
  double v[nbr_of_particles];
  double a[nbr_of_particles];

  v1 = malloc(sizeof(double)*avg_loops); v2 = malloc(sizeof(double)*avg_loops);
  v3 = malloc(sizeof(double)*avg_loops); v4 = malloc(sizeof(double)*avg_loops);
  x1 = malloc(sizeof(double)*avg_loops); x2 = malloc(sizeof(double)*avg_loops);
  x3 = malloc(sizeof(double)*avg_loops); x4 = malloc(sizeof(double)*avg_loops);

  for (int i=0, j=0; i<nbr_of_simulations; i++, j+=nbr_of_particles){
    x_pos[i] = j + x_pointer;
    v_particles[i] = j + v_pointer;
  }
  // First initialization of the position, velocities and acceleration
  for (int k=0; k<nbr_of_particles; k++){
          
    uni_var1 = gsl_rng_uniform(q);
    uni_var2 = gsl_rng_uniform(q);
    r_v  = sqrt( -2.0*log(uni_var1))*cos(M_PI*2*uni_var2);
    r_v = r_v * sigma;

    x_pos[0][k] = start_x;
    x[k] =  start_x;
    v_particles[0][k] = start_v;
    v[k] =  start_v;
    a[k] = get_acc(omega_0, x[k],m, eta,  r_v, v[k]);
  }
  // Loop to write a few tradjectories to disk
  for (int i=1; i<nbr_of_simulations; i++){

    for (int j=0; j<nbr_of_particles; j++){
      //gaussian_nbr(gauss_rv);

      uni_var1 = gsl_rng_uniform(q);
      uni_var2 = gsl_rng_uniform(q);
      uni_var3 = gsl_rng_uniform(q);
      uni_var4 = gsl_rng_uniform(q);

      // Using Box-Muller
      gauss_rv[0]  = sqrt( -2.0*log(uni_var1))*cos(M_PI*2*uni_var2);
      gauss_rv[1]  = sqrt( -2.0*log(uni_var1))*sin(M_PI*2*uni_var2);  
      r_v  = sqrt( -2.0*log(uni_var3))*cos(M_PI*2*uni_var4);
      r_v = r_v * sigma; // Change from unit variance to sigma² variance
      //v[j] = sqrt(c_0)*v[j] + v_th*sqrt(1.0-c_0)*gauss_rv[0];
      v[j] = 0.5*a[j]*timestep+sqrt(c_0)*v[j] + v_th*sqrt(1.0-c_0)*gauss_rv[0];
      x[j] = x[j] + v[j]*timestep;
      a[j] = get_acc(omega_0, x[j],m, eta,  r_v, v[j]);
      v[j] = sqrt(c_0)*0.5*a[j]*timestep + sqrt(c_0)*v[j] + v_th*sqrt(1.0-c_0)*gauss_rv[1];
      //v[j] = sqrt(c_0)*v[j] + v_th*sqrt(1.0-c_0)*gauss_rv[1];

      x_pos[i][j] = x[j];
      v_particles[i][j] = v[j];
    }
  }

  // Loop to calculate average over several trajectories
  for (int j=0; j<avg_loops; j++){

          tmp_v = start_v;
          tmp_x = start_x;

          uni_var1 = gsl_rng_uniform(q);
          uni_var2 = gsl_rng_uniform(q);
          r_v  = sqrt( -2.0*log(uni_var1))*cos(M_PI*2*uni_var2);
          r_v = r_v * sigma;

          tmp_a =  get_acc(omega_0, tmp_x, m, eta, r_v,tmp_v);
          
          v_mean[0] += tmp_v/avg_loops;
          v_variance[0] += tmp_v*tmp_v/avg_loops;
          
          x_mean[0] += tmp_x / avg_loops;
          x_variance[0] += tmp_x*tmp_x / avg_loops;

          for (int i=1; i<nbr_of_simulations; i++){

                  // Get four uniform random numbers [0,1]
                  
                  uni_var1 = gsl_rng_uniform(q);
                  uni_var2 = gsl_rng_uniform(q);
                  uni_var3 = gsl_rng_uniform(q);
                  uni_var4 = gsl_rng_uniform(q);

                  // Using Box-Muller
                  gauss_rv[0]  = sqrt( -2.0*log(uni_var1))*cos(M_PI*2*uni_var2);
                  gauss_rv[1]  = sqrt( -2.0*log(uni_var1))*sin(M_PI*2*uni_var2);
                  r_v  = sqrt( -2.0*log(uni_var3))*cos(M_PI*2*uni_var4);
                  r_v = r_v * sigma;

                  //tmp_v =  sqrt(c_0)*tmp_v + v_th*sqrt(1.0-c_0)*gauss_rv[0];
                  tmp_v = 0.5*tmp_a*timestep + sqrt(c_0)*tmp_v + v_th*sqrt(1.0-c_0)*gauss_rv[0];
                  tmp_x = tmp_x + tmp_v*timestep;
                  tmp_a =  get_acc(omega_0, tmp_x,m, eta,  r_v, tmp_v);
                  tmp_v = sqrt(c_0)*0.5*tmp_a*timestep + sqrt(c_0)*tmp_v + v_th*sqrt(1.0-c_0)*gauss_rv[1];
                  //tmp_v = sqrt(c_0)*tmp_v + v_th*sqrt(1.0-c_0)*gauss_rv[1];
                  
                  get_histogram_values(j, i,tmp_v, tmp_x);
                  
                  v_mean[i]+= tmp_v/avg_loops;
                  v_variance[i] += tmp_v*tmp_v/avg_loops;
                  x_mean[i] += tmp_x / avg_loops;
                  x_variance[i] += tmp_x*tmp_x / avg_loops;
                  
          }
  }
  for (int i=0; i<nbr_of_simulations; i++){
          v_variance[i] = sqrt(v_variance[i] - v_mean[i]*v_mean[i]); // std
          x_variance[i] = sqrt(x_variance[i] - x_mean[i]*x_mean[i]); // std
  }

  save_x_to_disk(x_pos, nbr_of_simulations, timestep, nbr_of_particles, x_mean, x_variance);

  save_v_to_disk(v_particles, nbr_of_simulations, timestep, nbr_of_particles, v_mean, v_variance);

  save_histogram_values_to_disk(avg_loops);

  return 0;
}



double get_acc(double omega_0, double x, double m, double eta, double r_v, double v){
  // F = -kx  = -m omega² x = m *a --> a = - omega² * x
  double a;
  double f;
  f = -m*omega_0*omega_0*x - m*eta*v + m*r_v;
  a = f / m;
  return a;
}


// Will generate guassian random nbr from central limit therorem.
void gaussian_nbr(double * gauss_rv){
  // generate uniform random number
  double uni_var1, uni_var2;
  uni_var1 = (double) rand() / (double) RAND_MAX;
  uni_var2 = (double) rand() / (double) RAND_MAX;
  // Using Box-Muller
  gauss_rv[0]  = sqrt( -2.0*log(uni_var1))*cos(M_PI*2*uni_var2);
  gauss_rv[1]  = sqrt( -2.0*log(uni_var1))*sin(M_PI*2*uni_var2);
}

void save_x_to_disk(double ** x_pos, int nbr_of_simulations,  double timestep, int nbr_of_particles, double * x_mean, double * x_variance){
        FILE * positions;
        double nm = 1E9; // save in nano meter
        double ms = 1E3; // save timestep ms
        positions = fopen("pos.dat", "w");
        double time_save;
        for (int i=0; i<nbr_of_simulations; i++){
                time_save = timestep*i*ms;
                fprintf(positions, "%e \t",time_save);
                for (int j=0; j<nbr_of_particles; j++){
                        fprintf(positions, "%e \t",x_pos[i][j]*nm);
                }
                fprintf(positions, "%e \t %e \t %e \t \n", x_mean[i]*nm, (x_mean[i]-x_variance[i])*nm, (x_mean[i]+x_variance[i])*nm);
        }
        fclose(positions);
}

void save_v_to_disk(double ** v_particles, int nbr_of_simulations,  double timestep, int nbr_of_particles, double *v_mean, double * v_variance){
        FILE * velocities;
        double ms = 1E3; // save in millisecond
        velocities = fopen("velocities.dat", "w");
        double time_save;
        for (int i=0; i<nbr_of_simulations; i++){
                time_save = timestep*i*ms;
                fprintf(velocities, "%e \t ",time_save);
                for (int j=0; j<nbr_of_particles; j++){
                        fprintf(velocities, "%e \t ",v_particles[i][j]*ms);
                }
                fprintf(velocities, "%e \t %e \t %e \t \n", v_mean[i]*ms, (v_mean[i]-v_variance[i])*ms, (v_mean[i]+v_variance[i])*ms);
        }
        fclose(velocities);
}

void get_histogram_values(int j, int i, double tmp_v, double tmp_x){
  int t_sample[4] = {1000, 2000, 5000, 20000};
  double ms = 1; //1E3;
  double nm = 1;//1E9;
  if (i == t_sample[0]) {
    v1[counter1] = tmp_v*ms;
    x1[counter1] = tmp_x*nm;
    counter1 +=1;
  }

  if (i == t_sample[1]) {
    v2[counter2] = tmp_v*ms;
    x2[counter2] = tmp_x*nm;
    counter2 +=1;
  }

  if (i == t_sample[2]) {
    v3[counter3] = tmp_v*ms;
    x3[counter3] = tmp_x*nm;
    counter3 +=1;
  }

  if (i == t_sample[3]) {
    v4[counter4] = tmp_v*ms;
    x4[counter4] = tmp_x*nm;
    counter4 +=1;
  }
}

void save_histogram_values_to_disk(int avg_loops){
  FILE * x_hist; FILE * v_hist;
  x_hist = fopen("x_hist.dat", "w");
  v_hist = fopen("v_hist.dat", "w");

  for (int i=0; i<avg_loops; i++){
    fprintf(x_hist, "%e \t %e \t %e \t %e\n", x1[i], x2[i], x3[i], x4[i]);
    fprintf(v_hist, "%f \t %f \t %f \t %f\n", v1[i], v2[i], v3[i], v4[i]);
  }
}
