#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define SIZE 100

void affiche(int tab[SIZE])
{
  for(int i=0 ; i<SIZE ; i++)
    printf("%d ", tab[i]);
  printf("\n");
}

int tableau[SIZE];
int resultat[SIZE] = {0};
int compteur[SIZE] = {0};
int size, rank;

int main(int argc,char **argv)
{
  srand(time(NULL));
    
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int tab2[SIZE/size];

  /* met des valeurs arbitraires dans le tableau : */
  for(int i=0 ; i<SIZE ; i++)
    tableau[i]=rand()%SIZE - 10;

  /* affiche le tableau initial (à supprimer pour de très grandes valeurs
     de SIZE) */
  affiche(tableau);
    
    



  /* compte le nombre de valeurs inférieures à chaque nombre du tableau
     dans le tableau compteur (complexité : SIZE^2).
     C'est ceci qu'il faut paralléliser en priorité. */
    
    
  MPI_Scatter(tableau, SIZE/size, MPI_INT, tab2, SIZE/size, MPI_INT, 0, MPI_COMM_WORLD);
    
  for(int i=0 ; i<SIZE ; i++)
  {
    for(int n=0 ; n<SIZE ; n++)
      if(tableau[n] < tableau[i])
        compteur[i]++;
  }


  MPI_Gather(tab2, SIZE/size, MPI_INT, tableau, SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    
  MPI_Finalize();
    
  /* crée le tableau résultat (complexité : SIZE).
     Inutile de paralléliser ceci. */
  for(int i=0 ; i<SIZE ; i++)
  {
    int d;
    /* les zéros c'est bon, le résultat est initialisé avec des 0. */
    if(tableau[i]==0)
      continue;
    /* cette boucle se place sur le premier élément qui n'est pas encore à
       la bonne valeur (pas 0) */
    for(d=0 ; ; d++)
      if(resultat[compteur[i]+d] != tableau[i])
        break;
    resultat[compteur[i]+d] = tableau[i];
  }


  /* affiche le tableau résultat (à supprimer pour de très grandes valeurs
     de SIZE) */
  affiche(resultat);

  return(0);
}
