all: knn_Parallel_MPI knn_Parallel_OpenMp

knn_Parallel_MPI: knn_Parallel.c
	mpic++ -o test knn_Parallel.c -fopenmp
	#First Program (Coppy paste and enter no of processors) "mpirun -np <Number of Processor> ./test"

knn_Parallel_OpenMp: knn_Parallel_OpenMP.c
	g++ knn_Parallel_OpenMP.c -fopenmp -o main
	#Second Program (Coppy paste) "./main"

Clean:
	rm -rf knn_Parallel_MPI 
	rm -rf knn_Parallel_OpenMp 


