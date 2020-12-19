//findout minimum from sub
//access individual component from each thread and print it on screen
//g++ test_Parallel_Update.c -fopenmp -o main

#include <iostream>
#include <fstream>
#include <sstream> 
#include <math.h> 
#include <bits/stdc++.h>
#include <chrono> 
#include <omp.h>

int	thread_count = 4;
using namespace std;
using namespace std::chrono; 

double a [150];
double b [150];
double c [150];
double d [150];

double X;
double Y;

double distance_A[90];
int testing = 99;//select any number between 46-50 (setosa) or 95-150 (Versicolor) 
				 //or 100-150(Verginica which is act as a random points)
int Index_array[90];

float Sequential_Time = 0; 
float Parallel_Time = 0; 
float Efficiency = 0; 
float Speed_up = 0;
float Serial_Fraction = 0; 

void Data_Sort(double a[], double b[], double c[], double d[]);
void EU_Distance(double b[], double d[], double X, double Y);
void Sorting(int Index_array[], double distance_A[]);
void* Odd_Even_Sort_2P(double *, int);

int main()
{
	cout << "K-NN Algorithm Using OpenMp - \n"<<endl;
	int k = 1; 
	int i, phase, index;
    double tmp;
	int n = 90;
	Data_Sort(a, b, c, d);
	double start_time = omp_get_wtime(); 
	EU_Distance(b, d, b[testing], d[testing]);

   for(phase = 0; phase < n; phase++)
   {
      if(phase % 2 == 0)	// Even phase, compare left
#        pragma omp parallel for num_threads(thread_count) \
            default(none) shared(distance_A, n, Index_array) private(i, index, tmp)
         for(i=1; i<n; i += 2)
         {
            if(distance_A[i-1] > distance_A[i])
            {
               tmp = distance_A[i-1];
               index = Index_array[i-1];
               distance_A[i-1] = distance_A[i];
               Index_array[i-1] = Index_array[i];
               distance_A[i] = tmp;
               Index_array[i] = index ;
            }
         }
      else			// Odd phase, compare right
#        pragma omp parallel for num_threads(thread_count) \
            default(none) shared(distance_A, Index_array, n) private(i, index, tmp)
         for(i=1; i<n-1; i += 2)
         {
            if(distance_A[i] > distance_A[i+1])
            {
               tmp = distance_A[i+1];
               index = Index_array[i+1];
               distance_A[i+1] = distance_A[i];
               Index_array[i+1] = Index_array[i];
               distance_A[i] = tmp;
               Index_array[i] = index ;
            }
         }
   }

   	double time = omp_get_wtime() - start_time;
    
	printf("---------------------------------------------------------------------------------------------------------\n");
	printf("Parallel Execution Time : %.2f Msec with %d Threads\n", 1000*time, thread_count+1);
	printf("Minimum Distance : %.2f Units\n", distance_A[0] );
	printf("Index of minimum Distance : %dth Out of 90 Test Samples\n", Index_array[0] );
	if (Index_array[0] >45)
	{
		printf("Test Sample (%.2f,%.2f) & Nearest Point to the Test Sample (%.2f,%.2f) from Versicolor\n", b[testing], d[testing],b[Index_array[0]], d[Index_array[0]]);
	}
	else
	{
		printf("Test Sample (%.2f,%.2f) & Nearest Point to the Test Sample (%.2f,%.2f) from Setosa\n", b[testing], d[testing],b[Index_array[0]], d[Index_array[0]]);
	}
	printf("---------------------------------------------------------------------------------------------------------\n");
	printf("Complte Sorted Array : \n");

    for (int i = 0; i < 90; ++i)
    {
    	cout << Index_array[i] << " : "<< distance_A[i] << "\t";
    }
}

void Sorting(int Index_array[],double distance_A[])
{

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

void* Odd_Even_Sort_2P(int *a, int n)
{
   int i, phase, tmp;

   for(phase = 0; phase < n; phase++)
   {
      if(phase % 2 == 0)	// Even phase, compare left
#        pragma omp parallel for num_threads(thread_count) \
            default(none) shared(a, n) private(i, tmp)
         for(i=1; i<n; i += 2)
         {
            if(a[i-1] > a[i])
            {
               tmp = a[i-1];
               a[i-1] = a[i];
               a[i] = tmp;
            }
         }
      else			// Odd phase, compare right
#        pragma omp parallel for num_threads(thread_count) \
            default(none) shared(a, n) private(i, tmp)
         for(i=1; i<n-1; i += 2)
         {
            if(a[i] > a[i+1])
            {
               tmp = a[i+1];
               a[i+1] = a[i];
               a[i] = tmp;
            }
         }
   }
   return NULL;
}  

void Data_Sort(double a[], double b[], double c[], double d[])
{
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