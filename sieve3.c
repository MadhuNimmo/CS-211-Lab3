/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;


   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

	low_value = 3 + 2 * (id*((n / 2) - 1) / p);
	high_value = 3 + 2 * (((id + 1)*(n / 2 - 1) / p) - 1);
	size = (high_value - low_value) / 2 + 1;
	cSize = 1000000;//Cache Size
	

	proc0_size = ((n / 2) - 1) / p;

	if ((3 + 2 * proc0_size) < (int)sqrt((double)n)) {
		if (!id) printf("Too many processes\n");
		MPI_Finalize();
		exit(1);
	}

	/* Allocate this process's share of the array. */

	marked = (char *)malloc(size);

	if (marked == NULL) {
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit(1);
	}


	low_0 = 3;
	high_0 = (int)sqrt((double)n);
	size_0 = (high_0 - low_0) / 2 + 1;
	marked_0 = (char *)malloc(size_0);
	index = 0;
	if (marked_0 == NULL) {
		printf("Cannot allocate enough memory for part2\n");
		MPI_Finalize();
		exit(1);
	}


	for (i = 0; i < size; i++) marked[i] = 0;
	for (i = 0; i < size_0; i++) marked_0[i] = 0;
	// Generating sieving primes
	index = 0;
	prime = 3;
	do {
		first_0 = (prime*prime - low_0) / 2;
		for (i = first_0; i<size_0; i += prime)
			marked_0[i] = 1;
		while (marked_0[++index]);
		prime = 2 * index + 3;

	} while (prime*prime <= n);
	loops = 0;
	marked = (char*)malloc(cSize);
	if (marked == NULL)
	{
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit(1);
	}

	if (!id)
		count = 1;
	else
	{
		count = 0;
	}
	prevLow = low_value;
	prevHigh = high_value;

	if (size%cSize == 0)
	{
		limit = size / cSize;
	}
	else
	{
		limit = size / cSize + 1;
	}
	do
	{
		low_value = ((prevLow)+(2 * loops*cSize));
		high_value = MIN((low_value + (cSize * 2 - 2)), prevHigh);
		nSize = (high_value - low_value) / 2 + 1;
		prime = 3;
		index = 0;
		for (i = 0; i<nSize; i++)
			marked[i] = 0;

		do
		{
			if (marked_0[index] == 0)
			{
				
				prime = 2 * index + 3;
				if (prime*prime > low_value)
					first = (prime*prime - low_value) / 2;
				else
				{
					if (!(low_value%prime))
						first = 0;
					else
					{
						if ((low_value%prime) % 2 != 0)
						{
							first = (prime - (low_value%prime)) / 2;
						}
						else
						{

							first = prime - (low_value%prime) / 2;
						}
					}
				}

				for (i = first; i<nSize; i += prime)
					marked[i] = 1;

			}
			index++;
		} while (index < (((high_0 - 3) / 2) + 1));
		int64_t j;
		for (j = 0; j<nSize; j++)
			if (!marked[j])
				count++;

		loops = loops + 1;
	} while (loops<limit); 

	if (p > 1) MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
		0, MPI_COMM_WORLD);

	
   /* Stop the timer */

   elapsed_time += MPI_Wtime();

   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

