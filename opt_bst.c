#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <values.h>
#include <stdlib.h>
#include <float.h>

#define MAXFLOAT FLT_MAX

void alloc_matrix(void ***, int, int, int ); 
void print_root(int **, int, int);
void print_matrix(int, int, void **);
void print_matrix_a(int, int, void **);

int main (int argc, char *argv[]) 
{
	FILE      *infileptr;
	float bestcost; /* Lowest cost subtree found so far */
	int bestroot; /* Root of the lowest cost subtree */
	int high; /* Highest key in subtree */
	int new_high;
	int i, j;
	int low; /* Lowest key in subtree */
	int new_low;
	int n; /* Number of keys */
	int r; /* Possible subtree root */
	float rcost; /* Cost of subtree rooted by r */
	int **root; /* Best subtree roots */
	float **cost; /* Best subtree costs */
	float *p; /* Probability of each key */

	//printf("Please give me the number of nodes:");
	//scanf("%d", &n);
	n = 500;
	p = (float *) malloc (n * sizeof(float));

	infileptr = fopen (argv[1], "r");
	fread(p, sizeof(float), n, infileptr);	

	
	/* find optimal binary search tree */
	alloc_matrix((void ***) &cost, n+1, n+1, sizeof(float));
	alloc_matrix((void ***) &root, n+1, n+1, sizeof(int));

	for (low = n; low >= 0; low--) 
	{
		cost[low][low] = 0.0;
		root[low][low] = low;

		for (high = low+1; high <= n; high++) 
		{
			bestcost = MAXFLOAT;
			/* distribute along the processors some of the calculations */
			if(high == low + 1)
			{
				bestcost = p[low];
				bestroot = root[low][low];
			}
			else
			{
				for (r = low; r < high; r++) 
				{
					rcost = cost[low][r] + cost[r + 1][high];
					for (j = low; j < high; j++) rcost += p[j];
					if (rcost < bestcost) 
					{
						bestcost = rcost;
						bestroot = r;
					}
				}
			}		
			cost[low][high] = bestcost;
			root[low][high] = bestroot;
		}
	}
	//print_root(root, 0, n-1);
	//print_matrix(n+1, n+1, (void **)root);
	//print_matrix_a(n+1, n+1, (void **)cost);
	return 0;
}

void print_root (int **root, int low, int high) 
{
	fprintf(stdout, "Root of the tree spanning %d-%d is %d\n",
	low, high, root[low][high+1]);
	fflush (stdout);
	if (low < root[low][high+1]-1)
	print_root (root, low, root[low][high+1]-1);
	if (root[low][high+1] < high-1)
	print_root (root, root[low][high+1]+1, high);
}
/* Allocate a two dimensional array with 'm' rows and
'n' columns, where each entry occupies 'size' bytes */
void alloc_matrix (void ***a, int m, int n, int size) 
{
	int i;
	void *storage;
	storage = (void *) malloc (m * n * size);
	*a = (void **) malloc (m * sizeof(void *));
	for (i = 0; i < m; i++) 
	{
		(*a)[i] = storage + i * n * size;
	}
}

void print_matrix(int m, int n, void **root)
{
	int i,j;
	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
			printf("%d   ", ((int **)root)[i][j]);
		printf("\n");
	}
		
}
void print_matrix_a(int m, int n, void **root)
{
	int i,j;
	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
			printf("%f     ", ((float **)root)[i][j]);
		printf("\n");
	}
		
}
//0.16 0.13 0.06 0.08 0.07 0.17 0.05 0.28
