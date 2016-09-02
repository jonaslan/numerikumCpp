#include<mpi.h>

int main (int argc, char** argv)
{
//	MPI program: no sharing of memory space: each process gets their own copy and there's a need of explicit sending the messages between each process
//	Initialisation of MPI execution environment
	MPI_Init(&argc,&argv);
	
	int size; //total number of processes: set at runtime by the user through mpirun wrapper
	int rank; //variable that identifies each process
	int sum = 0; // summation variable
	int id = 0; //label of the message passed
	MPI_Status status; //initialisation of status variable which contains information about the success or failure of a message transmission

//	Set size and rank of group associated to the communicator MPI_COMM_WORLD
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
//	Instructions for master's thread
	if(!rank)
	{
		MPI_Send(&sum, 1, MPI_INT, rank+1, id, MPI_COMM_WORLD); // send first message
		MPI_Recv(&sum, 1, MPI_INT, size-1, id, MPI_COMM_WORLD, &status); // receive message sent from the last thread (rank = size - 1); blocked until then
		std::cout << "Final result: " << sum << std::endl; // display final result
		std::cout << "Sum of arithmetic series: " << size * (size - 1) / 2 << std::endl; // check result
	}
 
//	Instructions for the last thread 
	else if (rank == size-1)
	{
		MPI_Recv(&sum, 1, MPI_INT, rank-1, id, MPI_COMM_WORLD, &status); // Receive message from previous thread
		sum += rank;
		MPI_Send(&sum, 1, MPI_INT, 0, id, MPI_COMM_WORLD); // Send back the processed message to master thread 
	}

//	Instructions for the other threads
	else
	{
		MPI_Recv(&sum, 1, MPI_INT, rank-1, id, MPI_COMM_WORLD, &status); // receive message from previous thread
		sum += rank;
		MPI_Send(&sum, 1, MPI_INT, rank+1, id, MPI_COMM_WORLD); // send message to next thread/process
	}

//	Finalise execution environment	
	MPI_Finalize();
}
