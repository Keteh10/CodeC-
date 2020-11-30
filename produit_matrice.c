#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>

/* matrix size */
#define N 2500

// Random values are [0, MAX_VAL]
#define MAX_VAL 5

// Nomber of checks
#define NBCHECKS 10
// acceptable error (in check)
#define ERROR   1.e-20

#define DIFFTEMPS(a,b) (((b).tv_sec - (a).tv_sec) + ((b).tv_usec - (a).tv_usec)/1000000.)

/* global to avoid stack overflow */
float a[N][N],b[N][N],c[N][N];

int main(int argc, char **argv)
{
  struct timeval tv_init, tv_begin, tv_end;
    gettimeofday( &tv_init, NULL);

  /***************************************************************************/
  // initialization
    
  srand((unsigned int)time(NULL));
    
//float a_colonne[N],
  float a_ligne[N];
//float b_colonne[N], b_ligne[N];
  float c_ligne[N];
  float sum;
  int size, rank;
  int i,j = 0;
    
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
    
  
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
    {
      b[i][j]=(double)rand()/(double)(RAND_MAX/MAX_VAL);
      a[i][j]=(double)rand()/(double)(RAND_MAX/MAX_VAL);
    }
    

  /***************************************************************************/
  // compute
  gettimeofday( &tv_begin, NULL);

                
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      c[i][j] = 0.;
    
  MPI_Scatter(a,N*N/size,MPI_FLOAT, a_ligne,N*N/size,MPI_FLOAT,0,MPI_COMM_WORLD);
     
  MPI_Bcast(b, N*N, MPI_INT, 0, MPI_COMM_WORLD);
    
  for(int i=0 ; i<N ; i++)
    for(int k=0 ; k<N ; k++)
      for(int j=0 ; j<N ; j++)
        c[i][j] += a[i][k] * b[k][j];


  gettimeofday( &tv_end, NULL);
    
    
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      sum = sum + a_ligne[j] * b[j][i]; }
        c_ligne[i] = sum;
        sum = 0;
        }

  MPI_Gather(c_ligne, N*N/size, MPI_INT, c, N*N, MPI_INT, 0, MPI_COMM_WORLD);


  MPI_Finalize();

  gettimeofday( &tv_end, NULL);

  /***************************************************************************/
  // check some values
  int check_ok = 1;
  for(int checks=0 ; checks<NBCHECKS ; checks++)
  {
    int i = rand()%N;
    int j = rand()%N;
    float val = 0.;
    for(int k=0 ; k<N ; k++)
      val += a[i][k] * b[k][j];
    if(fabs(val - c[i][j]) > ERROR)
    {
      fprintf(stderr, "BAD RESULTS !");
      fprintf(stderr, " (value[%d][%d] = %g should be %g)\n",
              i, j, c[i][j], val);
      check_ok = 0;
    }
  }
  if(check_ok)
    fprintf(stderr, "Ok results :)\n");

  /***************************************************************************/
  /* execution times */
  printf("Init : %lfs, Compute : %lfs\n",
         DIFFTEMPS(tv_init,tv_begin),
         DIFFTEMPS(tv_begin,tv_end));

  return( 0 );
}
