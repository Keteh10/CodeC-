#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 1


int ma_somme_globale(int *sendbuf, int *recvbuf);

int ma_somme_globale(int *sendbuf, int *recvbuf){
    
    int sum = 0;
    
    for (rank=0){
        sum = sendbuf;
    }
    
    for (int i=0; i<size; i++){
        MPI_Send(sendbuf, SIZE, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Recv(recvbuf, SIZE, MPI_INT, i ,0, MPI_COMM_WORLD);
        sum = sum + sendbuf;
    }
    
    for (int i=0; i<size; i++){
        recvbuf = sum;
        MPI_Send(recvbuf, SIZE, MPI_INT,0,i,MPI_COMM_WORLD);
    }
    return (sum);
}

int main(int argc, char **argv){
    
    int rank, size;
    int sum1;
    int tab[10]={1,2,3,4,5,6,7,8,9,10};
    int tab2[10];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    sum1 = ma_somme_globale(tab, tab2);
    MPI_Allreduce(tab, tab2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Finalize();
    
}
