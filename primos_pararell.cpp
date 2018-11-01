#include <iostream>
#include <fstream>
#include "omp.h"
#include "mpi.h"

using namespace std;

void imprimePrimos(int* vetor, int n) {
	cout << "-------------------------------" << endl;
	for (int i = 0; i <= n; i++) {
		if(vetor[i] != 1) {
			cout << i << ": " << vetor[i] << endl;
		}
	}	
	cout << "-------------------------------" << endl;
}

int main(int argc, char** argv) {
	int size, rank, rc;
	MPI_Init(&argc,&argv); 
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ofstream arquivo("THUZA");
	cout << "Sou o processo " << rank << " de " << size << endl;
	arquivo << "Sou o processo " << rank << " de " << size << endl;

	if(rank == 0) {

		int n = 33;
		int tamanhoVetor = n+1;

		int tam = (n+1) / size;
		int fatias[2];
		for(int i = 0; i < size; i++) {
			fatias[0] = i*tam;
			fatias[1] = fatias[0] + tam;
			if(i == n){
				fatias[1] = n;
			}
			MPI_Send(&fatias, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		// int tamanhoParte = tamanhoVetor / size;
		// int tamanhoPrimeiro = tamanhoParte + (tamanhoVetor/size);
		// MPI_Send(&n, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	int n;
	MPI_Status status;
	int fatias[2];
	MPI_Recv(&fatias, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	cout << "minha parte: " << fatias[0] << " ate " << fatias[1] << endl;

	// int vetor[tamanhoVetor];
	// vetor[0] = 1;
	// vetor[1] = 1;
	// for(int i = 2; i <= n; i++) {
	// 	vetor[i] = 0;
	// }

	// for (int i = 0; i <= n; ++i) {
	// 	cout << i << ": " << vetor[i] << endl;
	// }

	// for(int k = 2; k <= n; k++) {
	// 	// TESTAR SE O K TA MARCADO
	// 	#pragma omp parallel for
	// 	for (int i = k*k; i <= n; i++) {
	// 		#pragma omp critical
	// 		{
	// 			cout << "THREAD " << omp_get_thread_num() << " calculando " << i << " % " << k << endl;
	// 		}
	// 		if(i % k == 0) {
	// 			vetor[i] = 1;
	// 		}
	// 	}
	// }
	// // printa so os primos
	// imprimePrimos(vetor, n);
	MPI_Finalize(); 
	cout << endl;
	return 0;

}