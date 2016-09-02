#include<mpi.h>

int main (int argc, char** argv)
{

	MPI_Init(&argc,&argv);
	
	int size;
	int rank;
	int sum = 0;
	int id = 0;
	MPI_Status status;
	
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
	
	if(!rank)
	{
		MPI_Send(&sum, 1, MPI_INT, rank+1, id, MPI_COMM_WORLD);
		MPI_Recv(&sum, 1, MPI_INT, size-1, id, MPI_COMM_WORLD, &status);
		std::cout << sum << std::endl;
	} 
	else if (rank == size-1)
	{
		MPI_Recv(&sum, 1, MPI_INT, rank-1, id, MPI_COMM_WORLD, &status);
		sum += rank;
		MPI_Send(&sum, 1, MPI_INT, 0, id, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&sum, 1, MPI_INT, rank-1, id, MPI_COMM_WORLD, &status);
		sum += rank;
		MPI_Send(&sum, 1, MPI_INT, rank+1, id, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
}