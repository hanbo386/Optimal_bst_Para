#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <values.h>
#include <stdlib.h>
#include <float.h>
#include <mpi.h>
#include "MyMPI.h"

#define MAXFLOAT FLT_MAX

void alloc_matrix(void ***, int, int, int ); 
void print_root(int **, int, int);
void print_matrix(int, int, void **);
void print_matrix_a(int, int, void **);

int main (int argc, char *argv[]) 
{
	FILE      *infileptr;
	int proc_size_x, proc_size_y, proc_first_x, proc_first_y;
	int send_low, send_high, send_size;
	int recv_low, recv_high, recv_size;
	float bestcost; /* Lowest cost subtree found so far */
	int bestroot; /* Root of the lowest cost subtree */
	int high; /* Highest key in subtree */
	int new_high;
	int i, j, k, l;
	int low; /* Lowest key in subtree */
	int new_low;
	int n; /* Number of keys */
	int r; /* Possible subtree root */
	float rcost; /* Cost of subtree rooted by r */
	int **root; /* Best subtree roots */
	float **cost; /* Best subtree costs */
	float *p; /* Probability of each key */
	
	float min_cost;
	double time_start, time_end,time_total;
	int min_root;

	int *best_root_receive, *best_root_send;
	float *best_cost_receive, *best_cost_send;
	int *offset;
	int *size;
	
	typedef struct 
	{
		float val;
		float temp_val;
		int rank;
	}best_struct;
	
	//struct best_struct best_buffer;

	best_struct *best_send, *best_receive;
	MPI_Datatype best;
	MPI_Datatype old_types[2] = {MPI_FLOAT, MPI_INT};
	MPI_Aint indices[2], extent;

	int blocklens[2] = {2,1};
	
	int rank; 					/* variable that holds the id of the current process */
	int numprocs; 			/* variable that holds the total number os processors */
	n = 4000;
	
	/* Default MPI initialization */

	MPI_Status status;
	MPI_Init(&argc, &argv);
		
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	//if(!rank)
		//printf("%d \n", numprocs);
	indices[0] = 0;
	MPI_Type_extent(MPI_FLOAT, &extent);
	indices[1] = 2 * extent;

	MPI_Type_struct(2, blocklens, indices, old_types, &best);
	MPI_Type_commit(&best);	

	

	/*if(!rank)
	{
		
		printf("Please type in the number of keys:");
		scanf("%d", &n);
	}*/
	MPI_Barrier(MPI_COMM_WORLD);
	time_start = MPI_Wtime();
	/* broadcast to the other processors n's value */
	//MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	p = (float *) malloc (n * sizeof(float));
	if(!rank)
	{
		infileptr = fopen (argv[1], "r");
		fread(p, sizeof(float), n, infileptr);
	}
	/* broadcast to the other processors p's array value */
	MPI_Bcast (p, n, MPI_FLOAT, 0, MPI_COMM_WORLD);	

	
	/* find optimal binary search tree */
	alloc_matrix((void ***) &cost, n+1, n+1, sizeof(float));
	alloc_matrix((void ***) &root, n+1, n+1, sizeof(int));

	best_send = (best_struct *) malloc (n * n * sizeof(best_struct));
	best_receive = (best_struct *) malloc (n * n * sizeof(best_struct));
	//best_send = (int *) malloc (n * n * sizeof(int));
	/*best_root_send = (int *) malloc (n * n * sizeof(int));
	best_cost_send = (float *) malloc (n * n * sizeof(float));
	best_root_receive = (int *) malloc (n * n * sizeof(int));
	best_cost_receive = (float *) malloc (n * n * sizeof(float));*/
	
	int num = n+1;
	for(i = 1; i <= numprocs; i++)														//����numprocs������ʱ���򽫾����numprocs�ν��м���
	{
		
		/*����ÿ�����̻��ֵ������ݿ���ʼ�����Լ���ά��������*/
		proc_first_x = BLOCK_LOW(rank + i - 1, numprocs, num);
		proc_first_y = BLOCK_LOW(rank, numprocs, num);


		proc_size_x = BLOCK_SIZE(rank + i - 1, numprocs, num);
		proc_size_y = BLOCK_SIZE(rank, numprocs, num);


		//if(!rank)
		//	printf("%d, %d, %d, %d \n", proc_first_x, proc_first_y, proc_size_x, proc_size_y);

		/*�����������棬���㵱��Ҫ�������ݴ���ʱ����������Ҫ���ݵ����ݿ�Ĳ���*/
		send_low = BLOCK_LOW(rank, numprocs, num);
		send_high = BLOCK_HIGH(rank + i -1, numprocs, num) ;
		send_size = (send_high - send_low + 2)*(send_high - send_low + 1)/2;

		recv_low = BLOCK_LOW(rank + 1, numprocs, num);
		recv_high = BLOCK_HIGH(rank + i, numprocs, num);
		recv_size = (recv_high - recv_low + 2)*(recv_high - recv_low + 1)/2;

		//if(rank)
			//printf("%d, %d \n", send_low, send_high);

		/*���ڵ�һ�μ���ͺ���ļ��㻮�����ݷ�����ͬ���������ֶԴ�����һ�μ���ʱ���ֵ����ݿ�Ϊ��������״��ͬʱ��һ����Ҫ���жԽ������ݵĳ�ʼ��*/
		if(i == 1)															
		{
			for(low = proc_first_y + proc_size_y -1  ; low >= proc_first_y; low--)
			{
				root[low][low] = low;
				cost[low][low] = 0.0;
				
				for(high = low + 1; high < proc_first_y + proc_size_y; high++)
				{
					cost[n - low][n - high] = cost[n - low - 1][n - high] + p[low];
					bestcost = MAXFLOAT;
			
					new_low = root[low][high - 1];
					new_high = root[low + 1][high] + 1;
			
					if(high == low + 1)
					{
						bestcost = p[low];
						bestroot = root[low][low];
					}
					else
					{
						for (r = new_low; r < new_high; r++) 
						{
							rcost = cost[low][r] + cost[r + 1][high];
							rcost = rcost + cost[n - low][n - high];
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
		}
		
		/*���������ݻ���Ϊ����ʱ���еļ���*/
		else if(rank <= numprocs-i)																											
		{
			for(low = proc_first_y +proc_size_y - 1; low >= proc_first_y; low--)
			{	
				for(high = proc_first_x; high < proc_first_x + proc_size_x; high++)					//�˴����������ظ������պ���о����Ż�
				{
					//printf("fuck,%d,%d ", low, high);
					cost[n - low][n - high] = cost[n - low - 1][n - high] + p[low];
					bestcost = MAXFLOAT;
			
					new_low = root[low][high - 1];
					new_high = root[low + 1][high] + 1;
						
					if(high == low +1)
					{
						bestcost = p[low];
						bestroot = root[low][low];
					}
					else
					{
						for (r = new_low; r < new_high; r++) 
						{
							rcost = cost[low][r] + cost[r + 1][high];
							rcost = rcost + cost[n - low][n - high];
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
		}
		
		/*��0�Ž����⣬������̶���Ҫ����һ���̴��ݱ����̵ļ���������ÿ�ּ���֮����һ�ּ���ʱ���һ�����̽����������㡣*/
		if(rank && (numprocs-i >= rank ))									//�˴�Ϊ�ų�0�Ž����Ѿ�ÿ�ּ��㱻��̭���Ľ���														
		{
			l = 0;
			for(k = send_low; k <= send_high; k++)
			{
				for(j = k; j <= send_high; j++)
				{
					
					best_send[l].val = cost[k][j];
					best_send[l].temp_val = cost[n - k][n - j];
					best_send[l].rank = root[k][j];
					l++;
				}
			}
			/*l = 0;
			for(k = send_low; k <= send_high; k++)
			{
				for(j = send_low; j <= send_high; j++)
				{
					
					//best_root_send[l] = root[k][j];
					//best_cost_send[l] = root[k][j];
					best_send[l].temp_val = cost[n - k][n - j];
					best_send[l].rank = root[k][j];
					l++;
				}
			}*/
			MPI_Send(best_send, send_size, best, rank-1, 1, MPI_COMM_WORLD);
			//MPI_Send(best_root_send, send_size, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
			//MPI_Send(cost, n*n, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD);
			//printf("hi_1, %d \n", rank);
		}
		/*�����һ�������⣬������̶���Ҫ����һ���̽��������������ÿ�ּ���֮����һ�ּ���ʱ���һ�����̽����������㡣*/
		if((rank != numprocs-1)&&(numprocs-i > rank))
		{
			l = 0;
			MPI_Recv(best_receive, recv_size, best, rank+1, 1, MPI_COMM_WORLD, &status);
			//MPI_Recv(best_root_receive, recv_size, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
			//MPI_Recv(best_cost_receive, n*n, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);
			for(k = recv_low; k <= recv_high; k++)
			{
				for(j = k; j <= recv_high; j++)
				{
					//root[k][j] = best_root_receive[l];
					//cost[k][j] = best_cost_receive[l];
					cost[n - k][n - j] = best_receive[l].temp_val;
					cost[k][j] = best_receive[l].val;
					root[k][j] = best_receive[l].rank;
					l++;
				}
			}
		}
	}
	
	if(!rank)
	{
		time_end = MPI_Wtime();
		time_total = time_end - time_start;
		printf("%f \n", time_total);
		//print_root(root, 0, n-1);
		//print_matrix(n+1, n+1, (void **)root);
		//print_matrix_a(n+1, n+1, (void **)cost);
		//time_end = MPI_Wtime();
		//time_totle = time_end - time_start;
		//printf("\n%f \n", time_totle);
	}
	
	MPI_Finalize();
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
