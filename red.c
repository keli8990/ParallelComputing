#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>  
#include<time.h>
#define red 1;
#define blue 2;
#define white 0;


int finished;
int n_itrs = 0;
int redcount = 0;
int bluecount = 0;
int n, tile, convergence, MAX_ITRS;
int cover;
int max_bondary, min_bondary;


void init_grid(int ***grid, int grid_highth, int gridSize) {
	
	/*Init Row*/
	*grid = malloc(sizeof(int*) * grid_highth);
	/*Init Column*/
	int i;
	for (i = 0; i < grid_highth; i++) {
		(*grid)[i] = (int *)malloc(gridSize*sizeof(int));
	}
}






/*Fill grid with random red, blue, white*/
void fill_grid(int ***grid,int grid_highth, int gridSize) {
	/*printf("Highth:%d\n",grid_highth);*/
	struct timeval tv;  
    struct timezone tz = {0, 0};  
    gettimeofday(&tv, &tz);  
  
    /* tv.tv_usec milliseconds */  
    int seed = tv.tv_usec % 65536;  
    
	 
	float max = 1.0;
	int x, y;
	srand(seed);
	for (x = 0; x < grid_highth; x++)
	{
		for (y = 0; y < gridSize; y++)
		{
			float val = ((float)rand() / (float)(RAND_MAX)) * max;

			if (val <= 0.33)
			{
				(*grid)[x][y] = 0;
			}
			else if (val <= 0.66)
			{
				(*grid)[x][y] = 1;

			}
			else
			{
				(*grid)[x][y] = 2;
			}
		}

	}
	


}



int main(int argc, char **argv) {

	int myid, numprocs, source, destination;
	int myid_world, numprocs_world;

	finished = 0;
	int **grid; 	/* grid[row][col] */


	MPI_Status status;
	MPI_Init(&argc, &argv); /* Used to send the command line argumenys to all procs */
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs_world); /* Initializes the number of procs in the group specified by mpirun */
	MPI_Comm_rank(MPI_COMM_WORLD, &myid_world);	 /* Initialize the rank of this process in the group */

	
	n = atoi(argv[1]);
	tile = atoi(argv[2]);
	MAX_ITRS = atoi(argv[4]);
	convergence = atoi(argv[3]);





	if (numprocs == 1) {
		init_grid(&grid, n, n);
		fill_grid(&grid, n ,n);

		while (!finished && n_itrs < MAX_ITRS) {
			n_itrs++;

			/* red color movement */
			int i, j;
			for (i = 0; i < n; i++) {
				if (grid[i][0] == 1 && grid[i][1] == 0) {
					grid[i][0] = 4;
					grid[i][1] = 3;
				}
				for (j = 1; j < n; j++) {
					if (grid[i][j] == 1 && grid[i][(j + 1) % n] == 0) {
						grid[i][j] = 0;
						grid[i][(j + 1) % n] = 3;
					}
					else if (grid[i][j] == 3)
						grid[i][j] = 1;
				}
				if (grid[i][0] == 3)
					grid[i][0] = 1;
				else if (grid[i][0] == 4)
					grid[i][0] = 0;
			}

			/* blue color movement */
			for (j = 0; j < n; j++) {
				if (grid[0][j] == 2 && grid[1][j] == 0) {
					grid[0][j] = 4;
					grid[1][j] = 3;
				}

				for (i = 1; i < n; i++) {
					if (grid[i][j] == 2 && grid[(i + 1) % n][j] == 0) {
						grid[i][j] = 0;
						grid[(i + 1) % n][j] = 3;
					}
					else if (grid[i][j] == 3)
						grid[i][j] = 2;
				}
				if (grid[0][j] == 3)
					grid[0][j] = 2;
				else if (grid[0][j] == 4)
					grid[0][j] = 0;
			}

			/* count the number of red and blue in each tile and check if the computation can be terminated*/
			int tile_size = n / tile;
			int tile_total = tile_size * tile_size;
			int m, n;
			for (m = 0; m < tile_size; m++) {
				for (n = 0; n < tile_size; n++) {
					redcount = 0;
					bluecount = 0;
					int z, y;
					for (z = tile_size * m; z <tile_size*(m + 1); z++) {
						for (y = tile_size * n; y < tile_size*(n + 1); y++) {
							if (grid[z][y] == 1) {
								redcount++;
							}
							else if (grid[z][y] == 2)
							{
								bluecount++;
							}
						}
					}

					

					if ((((double)redcount) / tile_total) * 100 >= convergence) {
						cover = red;
						finished = 1;

						break;
					}
					else if ((((double)bluecount) / tile_total) * 100 >= convergence) {
						cover = blue;
						finished = 1;

						break;
					}
				}
				if (finished)
					break;
			}
			if (finished) {
				if (cover == 1) {
					printf("Red");
					printf("Block[%d][%d] \n", m, n);
					printf("Red takes  %d %% of the block \n", convergence);
					printf("iteration=%d \n", n_itrs);
					
					printf("Red count = %d\n", redcount);
					printf("Blue count = %d\n", bluecount);
					printf("Tiles size = %d\n", tile_total);
					printf("Red Percentage = %lf\n", (((double)redcount) / tile_total) * 100);
					printf("Blue Percentage = %lf\n\n\n", (((double)bluecount) / tile_total) * 100);
					break;
				}
				else if (cover == 2) {
					printf("Blue");
					printf("Block[%d][%d] \n", m, n);
					printf("Blue takes  %d %% of the tile \n", convergence);
					printf("iteration=%d \n", n_itrs);
					
					printf("Red count = %d\n", redcount);
					printf("Blue count = %d\n", bluecount);
					printf("Tiles size = %d\n", tile_total);
					printf("Red Percentage = %lf\n", (((double)redcount) / tile_total) * 100);
					printf("Blue Percentage = %lf\n\n\n", (((double)bluecount) / tile_total) * 100);
					break;
				}
			}


		}

		if (!finished) {
			printf("Not Converged");

		}
	}
	else {
		/*"tile" is the tile size input by user, the grid will be divided by "tile*tile" matrix*/
		int blocksize = n / tile;
		
		
		
		/*When the different processes have the same color, they belong to
		the same group, when the input process number exceed the require 
		(Tile)number, the maximum process number will be restricted in the 
		same group                                                     */
		int color = myid_world / tile;
		MPI_Comm row_comm;
		MPI_Comm_split(MPI_COMM_WORLD, color, myid_world, &row_comm);
		MPI_Comm_rank(row_comm, &myid);
		MPI_Comm_size(row_comm, &numprocs);
		if(color != 0){
			MPI_Barrier( row_comm);
		}
		else{
			printf("Numproc:%d\n",numprocs);
		destination = (myid + 1) % numprocs;
		source = (myid - 1 + numprocs) % numprocs;
		
		int cover;
		int i, j, k, l;
		
		
		
		
		int block_numpro = (int)(tile / numprocs);  
		int q = tile % numprocs;
		int *bufTop = malloc(sizeof(int) * n);
		int *bufBut = malloc(sizeof(int) * n);
		int *gather_Buf = malloc(sizeof(int) * numprocs);
		
		/*Calculating each process's max_boundary. when the number of processes could not
		be divided by the blocksize. The value of block_numpro = (int)(tile / numprocs)
		indicates that the regular number of lines of tiles for each process should handle. The value 
		q=tile % numprocs presents that how many process should be added one more line of tiles*/
		if (q != 0) {
				if (myid < q) {
					max_bondary = (myid + 1)*(block_numpro + 1)*blocksize - 1;
					min_bondary = myid*(block_numpro + 1)*blocksize;
				}
				else if (myid >= q) {
					min_bondary = myid * (block_numpro + 1)*blocksize - (myid - q)*blocksize;
					max_bondary = (myid + 1)*(block_numpro + 1)*blocksize - 1 - (myid - q + 1)*blocksize;
				}
			}
			else if (q == 0) {
				max_bondary = (myid + 1) * block_numpro * blocksize - 1;
				min_bondary = myid * block_numpro * blocksize;
			}
		
		int grid_Highth = max_bondary - min_bondary +1;
		init_grid(&grid, grid_Highth, n);
		
		fill_grid(&grid, grid_Highth, n);
		
		while (!finished && n_itrs < MAX_ITRS) {
			n_itrs++;

			
			
			/*   red movement    */
			for (i = 0; i < grid_Highth; i++) {
				if (grid[i][0] == 1 && grid[i][1] == 0) {
					grid[i][0] = 4;
					grid[i][1] = 3;
				}
				for (j = 1; j < n; j++) {
					if (grid[i][j] == 1 && (grid[i][(j + 1) % n] == 0)) {
						grid[i][j] = 0;
						grid[i][(j + 1) % n] = 3;
					}
					else if (grid[i][j] == 3) {
						grid[i][j] = 1;
					}
				}
				if (grid[i][0] == 3) {
					grid[i][0] = 1;
				}
				else if (grid[i][0] == 4) {
					grid[i][0] = 0;
				}

			}
			

			/*For each process, the value of Max Bondary pass the value to TopBuffer,
			Meanwhile, the value of Min Boundary will pass to bottomBuffer and save it.*/
			MPI_Sendrecv(grid[grid_Highth-1], n, MPI_INT, destination, 1, bufTop, n, MPI_INT, source, 1, row_comm, &status);

			MPI_Sendrecv(grid[0], n, MPI_INT, source, 1, bufBut, n, MPI_INT, destination, 1, row_comm, &status);
			
			
			/*The value stored in TopBuffer will be used to indicate whether the blue could move 
			in Min Boundary, if it could move, updating the Min Boundary's value        */
			int y = 0;
			for (y = 0; y < n; y++)
			{
				if (bufTop[y] == 2 && grid[0][y] == 0) {
					bufTop[y] = 4;
					grid[0][y] = 3;
				}
			}
			
			/*blue movement for each processor, restricted by Max and Min Boundary*/
			 /* blue color movement */
			for (j = 0; j < n; j++){
                if (grid[0][j] == 2 && grid[1][j] == 0){
                    grid[0][j] = 4;
                    grid[1][j] = 3;
                }
                for (i = 1 ; i < grid_Highth - 1; i++){
                    if (grid[i][j] == 2 && grid[(i+1)%n][j]==0){
                        grid[i][j] = 0;
                        grid[(i+1)%n][j] = 3;
                    }
                    else if (grid[i][j] == 3){
                        grid[i][j] = 2;
                    }
                }
                if (grid[0][j] == 3){
                    grid[0][j] = 2;
                }
                else if (grid[0][j] == 4){
                    grid[0][j] = 0;
                }
                
                
            }
			

			/*Comparing the value of BottomBuffer and the value in Max Boundary. Thus, determining whether 
			the Blue could move. If the blue could move, then, change the vcalue of Max Boundary*/
			int x;
			for (x = 0; x < n; x++)
			{
				if (bufBut[x] == 0 && grid[(max_bondary - min_bondary)][x] == 2) {
					bufBut[x] = 3;
					grid[(max_bondary - min_bondary)][x] = 4;
				}
			}
			
			/*Count Red and Blue in each peocess*/
			int m, n, s, p;

			for (m = 0; m < (grid_Highth /blocksize); m++) {
				for (n = 0; n < blocksize; n++) {
					redcount = 0;

					bluecount = 0;
					for (s = m*blocksize; s < (m + 1)*blocksize; s++) {
						for (p = n*blocksize; p < (n + 1)*blocksize; p++) {
							if (grid[s][p] == 1) {
								redcount++;
							}
							else if (grid[s][p] == 2) {
								bluecount++;
							}
						}
					}

					int tile_size = blocksize*blocksize;

					
					if ((((double)redcount) / tile_size) * 100 >= convergence) {

						cover = 1;
						finished = 1;
						
						break;
					}
					else if ((((double)redcount) / tile_size) * 100 >= convergence) {
						
						cover = 2;
						finished = 1;
						break;
					}
				}
				
				if (finished) {
					if (cover == 1) {
						printf("Block[%d][%d] \n", m, n);
						printf("Red count = %d\n", redcount);
						printf("Blue count = %d\n", bluecount);
						printf("Block size = %d\n",(blocksize * blocksize));
						printf("Red takes more than %d %% of the tile \n", convergence);
						printf("Proc %d already stopped, iteration: %d.\n", myid, n_itrs);
						printf("Red Percentage = %lf\n", (((double)redcount) / (blocksize * blocksize)) * 100);
						printf("Blue Percentage = %lf\n\n\n", (((double)bluecount) / (blocksize * blocksize)) * 100);
						break;
					}
					else if (cover == 2) {
						printf("Block[%d][%d] \n", m, n);

						printf("Red count = %d\n", redcount);
						printf("Blue count = %d\n", bluecount);
						printf("Block size = %d\n",(blocksize * blocksize));
						printf("Blue takes more than %d %% of the tile \n", convergence);
						printf("Proc %d already stopped, iteration: %d.\n", myid, n_itrs);
						printf("Red Percentage = %lf\n", (((double)redcount) / (blocksize * blocksize)) * 100);
						printf("Blue Percentage = %lf\n\n\n", (((double)bluecount) / (blocksize * blocksize)) * 100);
						break;
					}
				}
				
			}



			/*The MPI_Reduce could get all the process's finished's value, and do "LOR" 
			computation,Then, save the value in Allfinished. The MPI_Borcast could help
			to broadcast to the other processes*/
			
			int allFinished;
			MPI_Reduce(&finished, &allFinished, 1, MPI_INT, MPI_LOR, 0,row_comm);
			finished = allFinished;
			MPI_Bcast(&finished, 1, MPI_INT, 0, row_comm);


		}
		if (!finished)
		{
			printf("Iteration number has execced the max Iteration number\n", n_itrs);
		}

		free(grid);

	}

	}

	MPI_Finalize();
	return 0;
}
