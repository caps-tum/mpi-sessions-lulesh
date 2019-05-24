
#include <mpi.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <cstdlib>

void InitMeshDecomp(int numRanks, int myRank,
                    int *col, int *row, int *plane, int *side)
{
   int testProcs;
   int dx, dy, dz;
   int myDom;
   
   // Assume cube processor layout for now 
   testProcs = int(cbrt(double(numRanks))+0.5) ;
   if (testProcs*testProcs*testProcs != numRanks) {
      printf("Num processors must be a cube of an integer (1, 8, 27, ...)\n") ;
#if USE_MPI  
printf("MPI_ABORT called in %s",  __PRETTY_FUNCTION__);    
      MPI_Abort(mpi.current_comm, -1) ;
#else
      exit(-1);
#endif
   }
   if (sizeof(double) != 4 && sizeof(double) != 8) {
      printf("MPI operations only support float and double right now...\n");
#if USE_MPI      
printf("MPI_ABORT called in %s",  __PRETTY_FUNCTION__);
      MPI_Abort(mpi.current_comm, -1) ;
#else
      exit(-1);
#endif
   }

   dx = testProcs ;
   dy = testProcs ;
   dz = testProcs ;

   // temporary test
   if (dx*dy*dz != numRanks) {
      printf("error -- must have as many domains as procs\n") ;
#if USE_MPI      
printf("MPI_ABORT called in %s",  __PRETTY_FUNCTION__);
      MPI_Abort(mpi.current_comm, -1) ;
#else
      exit(-1);
#endif
   }
   int remainder = dx*dy*dz % numRanks ;
   //printf("Remainder %d", remainder);
   if (myRank < remainder) {
      myDom = myRank*( 1+ (dx*dy*dz / numRanks)) ;
   }
   else {
      myDom = remainder*( 1+ (dx*dy*dz / numRanks)) +
         (myRank - remainder)*(dx*dy*dz/numRanks) ;
   }
     // printf("myDom %d", myDom);


   *col = myDom % dx ;
   *row = (myDom / dx) % dy ;
   *plane = myDom / (dx*dy) ;
   *side = testProcs;

   return;
}


void pack(int numRanks,
                    int myRank,
                              std::vector<int>& input,
                              std::vector<int>& output
                             )
{
  
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    //Laik_Space* space = laik_partitioning_get_space(p);
    //const Laik_Slice* slice = laik_space_asslice(space);
    //int edgeElems= (int) (cbrt( (domain.sizeX()+1) / numRanks ) + 0.1 );

    int Nx=3,Ny=3,Nz=3;

    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*(Nx-1)+1;
    int Ly = Ry*(Ny-1)+1;
    int Pxy= Lx*Ly;
    //int Pyz= Ly*Lz;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this
    // task
    unsigned long from, to, n=0;
    int r=0;
    int nx=0;
    int tag=0;
    int d=1;

    printf("Rank: %d, Rx, Ry, Rz = %d, Lz = %d, Ly= %d, Pxy= %d\n", myRank, Rx, Lx, Ly, Pxy);

    for (int rz = 0; rz < Rz; rz++)
    {
        for (int ry = 0; ry < Ry; ry++)
        {
            for (int rx = 0; rx < Rx; rx++)
            {
                r = rx + ry*Rx + rz*Rx*Ry; // (xyz)
                if(r!=myRank) continue;

                //r = rx + rz*Rx + ry*Rx*Rz; // (xzy)
                //r = ry + rz*Ry + rx*Rz*Ry; // (yzx)
                //r = rz + ry*Rz + rx*Rz*Ry; // (zyx)
                //r = rz + rx*Rz + ry*Rx*Rz; // (zxy)

                // loop over y and z to create the slices in the
                // partitioning
                for (int ny = 0 ; ny < Ny; ny++)
                {
                    for (int nz = 0 ; nz < Nz; nz++)
                    {
                        nx=0;
                        // unique tags
                        tag = nx + Lx*ny + Pxy*nz +  rx*(Nx-1) + ry*Lx*(Ny-1) + rz*Pxy*(Nz-1) + Ny*100;
                        nx=0;
                        from = nx + Lx*ny + Pxy*nz +  rx*(Nx-1) + ry*Lx*(Ny-1) + rz*Pxy*(Nz-1);
                        nx=Nx;
                        to = nx + Lx*ny + Pxy*nz +  rx*(Nx-1) + ry*Lx*(Ny-1) + rz*Pxy*(Nz-1);
                        nx=0;
                        n = nx + ny*Nx + nz*Nx*Ny;
                        #pragma ivdep
                        for(int i=0; i<(to-from); i++){
                             output[i+from] = input[i+n];                
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[]){

        int size;
        int rank;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        int col, row, plane, side;
        InitMeshDecomp(size, rank, &col, &row, &plane, &side);

        std::vector<int> input(3*3*3, rank);
        std::vector<int> output(5*5*5,0);
        
        pack(size, rank, input, output);

        if(rank==0){
        for (int i=0; i<output.size(); i++){
            printf("OutVector: %d\n", output[i]);
        }
        }
        for (int i=0; i<input.size(); i++){
            printf("InputVector: %d\n", input[i]);
        }
        std::vector<int> soutput(5*5*5,0);

        MPI_Allreduce(output.data(), soutput.data(), 5*5*5, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(rank==0){
        for (int i=0; i<soutput.size(); i++){
            printf("ReducedVector: %d\n", soutput[i]);
        }
        }

        MPI_Finalize();
    
    }