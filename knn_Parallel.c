//findout minimum from sub
//access individual component from each thread and print it on screen
// mpic++ -g -Wall -o test test9.c -fopenmp  
//mpirun -np 5 ./test


/*-----------------------------------------------------------------
 * File:    knn_Parallel.c
 * 
 * Purpose: Knn Parallel Algorithm using using OpenMP and MPI
 *
 * Compile: mpic++ -g -Wall -o test knn_Parallel.c -fopenmp 
 *   
 * Run:     mpirun -np <Number of Processor> ./test
 *
 * Input:   2 Parameters (sepal width and petal width) 
 * Output:  Calculate Euclidian Distance using OpenMP and then pass 
 			it through MPI parallel loop to calculate shortest 
 			distance among all data samples. Findout the index of 
 			that data sample(which is not paralleliz yet)   
 			Performance Analysis Parameters: Speed_up, Efficiency, 
 			Serial_Fraction
 ------------------------------------------------------------------*/
#include <iostream>
#include <sstream> 
#include <math.h> 
#include <bits/stdc++.h>
#include <chrono> 
#include <omp.h>
#include <mpi.h> 
#include <algorithm>
 #define n 90  //Data Samples

using namespace std;
using namespace std::chrono; 

double a [150]; //Empty array to store Iris data set
double b [150];
double c [150];
double d [150];

double X;
double Y;

double distance_A[90];
int testing = 49;// For testing purpose : select any number between 46-50 (setosa) 
				 //or 95-150 (Versicolor) or 100-150(Verginica which is act as a random points)

int Index_array[90];
double a2[1000]; 
  
void Data_Sample(double a[], double b[], double c[], double d[]);
void EU_Distance(double b[], double d[], double X, double Y);
void Sorting(int Index_array[], double distance_A[]);
void parallel_loop();

int main(int argc, char* argv[]) 
{
int pid, np, // np - no. of processes and pid - process id 
        elements_per_process, 
        n_elements_recieved; 
    double elapsed_time; 

    MPI_Status status; // Creation of parallel processes 
    MPI_Init(&argc, &argv); // find out process ID and how many processes were started
     
    MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
    MPI_Comm_size(MPI_COMM_WORLD, &np); //Assign rank and processes

	int k = 1; //only 1 neighbour 

	int thread_num = 9;//Thread count for OpenMP

	Data_Sample(a, b, c, d); //Sort Data

 	omp_set_num_threads(thread_num); //OpenMP Parameter

    MPI_Barrier(MPI_COMM_WORLD); //Start all Processes at the same time

    elapsed_time = -MPI_Wtime();  //MPI Time Calculation

	EU_Distance(b, d, b[testing], d[testing]); //OpenMP as a parallel Processing Parameter
    
//	Sorting(Index_array, distance_A); //For Serial Testing Purpose

	// master process 
        if (pid == 0) { 
        		cout << "K-NN Algorithm Using OpenMp and MPI - "<<endl;
				cout << "-------------------------------------------------------------------------------------------"<<endl;
	
            int index, i; 
            elements_per_process = n / np; 
      
            // check if more than 1 processes are run 
            if (np > 1) { 
                // distributes the portion of array 
                // to child processes to calculate 
                // Get Minimum out of array
                for (i = 1; i < np - 1; i++) { 
                    index = i * elements_per_process; 
                    MPI_Send(&elements_per_process, 
                             1, MPI_INT, i, 0, 
                             MPI_COMM_WORLD); 
                    MPI_Send(&distance_A[index], 
                             elements_per_process, 
                             MPI_DOUBLE, i, 0, 
                             MPI_COMM_WORLD); 
                } 
      
                index = i * elements_per_process; 

                int elements_left = n - index; 
      
                MPI_Send(&elements_left, 
                         1, MPI_INT, 
                         i, 0, 
                         MPI_COMM_WORLD); 
                MPI_Send(&distance_A[index], 
                         elements_left, 
                         MPI_DOUBLE, i, 0, 
                         MPI_COMM_WORLD); 
            } 
            // master process get its minimum from its sub array 
            double sum[elements_per_process]; 
            for (i = 0; i < elements_per_process; i++) 
                {
                	sum[i] = a[i];
                } 
            double small = *min_element(sum, sum + elements_per_process);
            // collects all minimum values from other processes 
            double tmp; 
            
            int new_count = np;
            double tmp1[new_count];

            for (i = 1; i < np; i++) 
            { 
                MPI_Recv(&tmp, 1, MPI_DOUBLE, 
                         MPI_ANY_SOURCE, 0, 
                         MPI_COMM_WORLD, 
                         &status); 
                int sender = status.MPI_SOURCE; 
                tmp1[i] = tmp;
            } 
            tmp1[0] = small;

            sort(tmp1,tmp1 + new_count);
                elapsed_time += MPI_Wtime();
                int distance_index = 0;
                double parallel_Time = 1000*elapsed_time;

                cout << "Parallel Execution Time : "<<parallel_Time << " Msec with " << np << endl;
               auto itr = find(distance_A, distance_A+90, tmp1[0]);
               distance_index = distance(distance_A, itr) + 5; //adjusting distance by saperating 
               //5 values for testing purpose 
                cout << "Minimum Distance : " << distance_A[distance_index]<< "Units"<< endl;
               cout << "Index of minimum distance : "<<distance_index <<" Out of 90 Test Sample"<<endl;
              if (distance_index > 45)
              {
              cout << "Test Sample ("<< b[testing] << ","<< d[testing] <<") & Nearest Point to the Test Sample (" << b[distance_index+5] << ","<< d[distance_index+5] <<") from Versicolor"<<endl; 
               }
              else
              {
              cout << "Test Sample ("<< b[testing] << ","<< d[testing] <<") & Nearest Point to the Test Sample (" << b[distance_index] << ","<< d[distance_index] <<") from Setosa"<<endl; 
               }	

				cout << "-----------------------------------------------------------------------------------------"<<endl;
 
        } 
    // slave processes 
    else { 
        MPI_Recv(&n_elements_recieved, 
                 1, MPI_INT, 0, 0, 
                 MPI_COMM_WORLD, 
                 &status); 
  
        // stores the received array segment 
        // in local array a2 
        MPI_Recv(&a2, n_elements_recieved, 
                 MPI_DOUBLE, 0, 0, 
                 MPI_COMM_WORLD, 
                 &status); 
  
        // Findout Minimum
        double temp[n_elements_recieved];
        double test; 
        for (int i = 0; i < n_elements_recieved; i++) 
        {
            temp[i] = a2[i];
        }
        test = *min_element(a2, a2 + n_elements_recieved);
        MPI_Send(&test, 1, MPI_DOUBLE, 
                 0, 0, MPI_COMM_WORLD); 
    } 
    // cleans up all MPI state before exit of process 
    MPI_Finalize(); 
}

void Sorting(int Index_array[],double distance_A[])
{
    //Serial Sort function with indexing for testing purpose
	pair<double, int> pairt[90];
	for (int i = 0; i < 90; i++)
	{
		pairt[i].first = distance_A[i];
		pairt[i].second = Index_array[i];
	}
	sort(pairt, pairt + 90);
	for (int l = 0; l < 90; l++)
	{
		distance_A[l] = pairt[l].first;
		Index_array[l] = pairt[l].second;
	}
}
void EU_Distance(double b[], double d[], double X, double Y)
{
	//Calculate EU Distance between test point and all sample points 
	//Also included Index
	double temporary;

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < 45; i++)
	{
		temporary = pow((b[i] - X),2) + pow((d[i] - Y),2) ;
		distance_A[i] = sqrt(temporary);
		temporary = 0;
		temporary = pow((b[50+i] - X),2) + pow((d[50+i] - Y),2) ;
		distance_A[45+i] = sqrt(temporary);
		temporary = 0;
		Index_array[i] = i; 
		Index_array[45+i] = 50+i;
	}
}
void Data_Sample(double a[], double b[], double c[], double d[])
{
	//Convert data into array from .txt file
	//Save in 4 arrays
	string myText;
	string newstr;
	string testing;
	char temp; 
	int counter = 0;
	ifstream MyReadFile("iris.txt");

	while (getline (MyReadFile, myText)) 
	{
		for(int i = 0; i<myText.size(); i++)
			{
			 if(myText[i] == ',')
				{
				 	temp = myText[i-3];
				 	newstr += temp;
				 	temp = myText[i-2];
				 	newstr += temp;
				 	temp = myText[i-1];
				 	newstr += temp;
				}
		}
	}

	MyReadFile.close();

	for (int j = 0; j < newstr.length(); j+=12)
		{
			testing += newstr[0+j];
			testing += newstr[1+j];
			testing += newstr[2+j];
			a[counter] = stof(testing);
			testing.clear();
			testing += newstr[3+j];
			testing += newstr[4+j];
			testing += newstr[5+j];
			b[counter] = stof(testing);
			testing.clear();
			testing += newstr[6+j];
			testing += newstr[7+j];
			testing += newstr[8+j];
			c[counter] = stof(testing);
			testing.clear();
			testing += newstr[9+j];
			testing += newstr[10+j];
			testing += newstr[11+j];
			d[counter] = stof(testing);
			testing.clear();
			counter = counter + 1;
		}		
}