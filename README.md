# Parallel_Programing_KNN

# Introduction
K-nearest neighbors algorithm is unsupervised learning algorithm which is used
to allocate input data into nearest cluster. This algorithm also called slow
algorithm because of its slow scanning nature. It can also used for regression or
classification purpose

# Algorithm and libraries
For this project I used Iris data set. Each row of the table represents an iris
flower, including its species and dimensions of its botanical parts, sepal and
petal, in centimeters.
Note - 1) I prefer to use 50 data points of setosa and 50 data points from
1versicolor from Iris data-set. Out of 50 samples i used 45 as a training data
and 5 as a testing data. User can select respective test sample by accessing
knn Parallel.c code and change parameter ”testing”. So in total 90 training
data and 10 testing data.
2) I prefer to use K = 1.
3) Euclidean Norm is used

# Libraries :
1) stdio.h: File which contains C declaration and Macro definition to be
shared between several files. stdio.h means standard input/output function
which contains printf(), scanf() functions.
2) limits.h: This header defines constants with the limits of fundamental
integral types for the specific system and compiler implementation used.
3) omp.h: Library for parallel programming in the SMP (symmetric multi-
processors, or shared-memory processors) model.
4) time.h This header file contains definitions of functions to get and ma-
nipulate date and time information.
5) stdlib.h This header defines several general purpose functions.
6) omp.h: Library for parallel programming in the SMP (symmetric multi-
processors, or shared-memory processors) model.
7) iostream: It stands for standard input-output stream. This header file
contains definitions to objects like cin, cout, cerr etc.
8) sstream: Its a stringstream class in C++ is a Stream Class to Operate
on strings.
9) bits/stdc++.h: Its a header file. This file includes all standard library.
10) math.h: eclares a set of functions to compute common mathematical
operations and transformations.
11) chrono: It deal with time. This is done using three concepts: Dura-
tion’s, Time points, Clocks.

# Serial Approach :
1) Data being loaded from .text file. At the beginning of the program it initial-
ize 4 different arrays with the size of 150 each. I used ”,” separation approach
and extract 4 parameters into arrays.
2) As per the note mentioned I took 90 samples with testing data and find-out
euclidean norm with correct indexing of training data.
3) Further data sorted-out along with index in ascending order. This sorting
include data sort along with indexing.
4) My solution uses 1-NN, the first data point will have minimum distance from
testing data sample to nearest neighbor.

# Parallel Approach :
There are two different parallel program I made. The first program is based
on OpenMP and second program uses OpenMP and MPI. I followed Fosters
approach in my second solution. I have attached both test results in the result
section. For now i am only explaining second solution.

# Foster Approach :
## Partition:
The program consist of 3 different type of threads. First type of threads is
used to calculate EU distance with respect to testing data point. I assign 10
threads to calculate EU distance and store everything in array EU Array[]. Sec-
ond type of thread act as a Master thread and slice EU Array[] into multiple
pieces. Number of pieces are equal to assign threads - 1. It send array slices to
all Slave threads.Third type of threads are nothing but Slave threads which are
used to find out minimum number out of its array.
## Communication:
Overall Program has communication link between EU array[], Master thread[]
and Slave thread[]. Data of EU array[] is accesses by Master thread[] then Mas-
ter thread[] will establish communication with other Slave threads[].
## Agglomeration:
The result of EU Array[] which is nothing but thread number 0 is communicate
with master thread Master thread[] is grouped with each slave thread
will execute its own command and at the end send its index to Master thread[].
Its a fix procedure entire program will follow this.
## Mapping:
EU distance calculation is divide into multiple threads which is assign in pro-
gram itself. In second and third type of task’s thread controlled by user. array
points are equally divide assign to Slave thread[]. If there is any extra points
left then master process will compute its minimum value.
