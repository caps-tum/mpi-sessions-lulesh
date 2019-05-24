#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

#if USE_MPI
#include <mpi.h>
#endif
#include "lulesh.h"

/* Helper function for converting strings to ints, with error checking */
template<typename IntT>
int StrToInt(const char *token, IntT *retVal)
{
   const char *c ;
   char *endptr ;
   const int decimal_base = 10 ;

   if (token == NULL)
      return 0 ;
   
   c = token ;
   *retVal = strtol(c, &endptr, decimal_base) ;
   if((endptr != c) && ((*endptr == ' ') || (*endptr == '\0')))
      return 1 ;
   else
      return 0 ;
}



void mapLocaltoGlobalNodes(int numRanks,
                              int myRank,
                              Domain& locDom,
                              std::vector<double>& m_fx,
                              std::vector<double>& m_fy,
                              std::vector<double>& m_fz,
                              std::vector<double>& m_x,
                              std::vector<double>& m_y,
                              std::vector<double>& m_z,
                              std::vector<double>& m_xd,
                              std::vector<double>& m_yd,
                              std::vector<double>& m_zd,
                              std::vector<double>& m_xdd,
                              std::vector<double>& m_ydd,
                              std::vector<double>& m_zdd,
                              std::vector<double>& m_nodalMass)
{
   int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    //Laik_Space* space = laik_partitioning_get_space(p);
    //const Laik_Slice* slice = laik_space_asslice(space);
    //int edgeElems= (int) (cbrt( (domain.sizeX()+1) / numRanks ) + 0.1 );

    int Nx=locDom.sizeX()+1,Ny=locDom.sizeX()+1,Nz=locDom.sizeX()+1;

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


    //printf("Rank: %d, Rx, Ry, Rz = %d, Lz = %d, Ly= %d, Pxy= %d\n", myRank, Rx, Lx, Ly, Pxy);

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
                           m_fx[i+from] = locDom.m_fx[i+n];
                           m_fy[i+from] = locDom.m_fy[i+n];
                           m_fz[i+from] = locDom.m_fz[i+n];
                           m_nodalMass[i+from] = locDom.m_nodalMass[i+n];
                           m_x[i+from] = locDom.m_x[i+n];
                           m_y[i+from] = locDom.m_y[i+n];
                           m_z[i+from] = locDom.m_z[i+n];
                           m_xd[i+from] = locDom.m_xd[i+n];
                           m_yd[i+from] = locDom.m_yd[i+n];
                           m_zd[i+from] = locDom.m_zd[i+n];
                           m_xdd[i+from] = locDom.m_xdd[i+n];
                           m_ydd[i+from] = locDom.m_ydd[i+n];
                           m_zdd[i+from] = locDom.m_zdd[i+n];             
                        }
                    }
                }
            }
        }
    }
}

void mapLocaltoGlobalElements(int numRanks,
                              int myRank,
                              Domain& locDom,
                              std::vector<double>& m_delv_xi,
                              std::vector<double>& m_delv_eta,
                              std::vector<double>& m_delv_zeta,
                              std::vector<double>& m_dxx,
                              std::vector<double>& m_dyy,
                              std::vector<double>& m_dzz,
                              std::vector<double>& m_delx_xi,
                              std::vector<double>& m_delx_eta,
                              std::vector<double>& m_delx_zeta,
                              std::vector<double>& m_e,
                              std::vector<double>& m_p,
                              std::vector<double>& m_q,
                              std::vector<double>& m_ql,
                              std::vector<double>& m_qq,
                              std::vector<double>& m_v,
                              std::vector<double>& m_volo,
                              std::vector<double>& m_delv,
                              std::vector<double>& m_vdov,
                              std::vector<double>& m_arealg,
                              std::vector<double>& m_ss,
                              std::vector<double>& m_elemMass
                             )
{
  
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    //Laik_Space* space = laik_partitioning_get_space(p);
    //const Laik_Slice* slice = laik_space_asslice(space);
    //int edgeElems= (int) (cbrt( (domain.sizeX()+1) / numRanks ) + 0.1 );

    int Nx=locDom.sizeX(),Ny=locDom.sizeX(),Nz=locDom.sizeX();
    int Nxh=locDom.sizeX()+1, Nyh=locDom.sizeX()+1, Nzh=locDom.sizeX()+1;

    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*Nx;
    int Ly = Ry*Ny;
    //int Lz = Rz*Nz;
    int Pxy= Lx*Ly;
    //int Pxz= Lx*Lz;
    //int Pyz= Ly*Lz;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this
    // task
    unsigned long from, to, n;
    int r=0;
    int nx=0;
    int tag=0;
    for (int rz = 0; rz < Rz; rz++)
    {
        for (int ry = 0; ry < Ry; ry++)
        {
            for (int rx = 0; rx < Rx; rx++)
            {
                r = rx + ry*Rx + rz*Rx*Ry; // task number
                if(r!=myRank) continue;

                // loop over z and x  to create the slices in the
                // partitioning
                for (int ny = 0; ny < Ny; ny++)
                {
                    for (int nz = 0; nz < Nz; nz++)
                    {   
                        // tag = global index where nx = 0 + safety shift = Ny+10
                        nx=0;
                        tag = nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz + Ny*10;
                        nx=0;
                        from = nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx=Nx;
                        to = nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx=0;
                        n = nx + ny*Nx + nz*Nx*Ny;
                        #pragma ivdep
                        for(int i=0; i<(to-from); i++){
                             m_delv_xi[i+from] = locDom.m_delv_xi[i+n];
                             m_delv_eta[i+from] = locDom.m_delv_eta[i+n];
                             m_delv_zeta[i+from] = locDom.m_delv_zeta[i+n];
                             m_dxx[i+from] = locDom.m_dxx[i+n];
                             m_dyy[i+from] = locDom.m_dyy[i+n];
                             m_dzz[i+from] = locDom.m_dzz[i+n];
                             m_delx_xi[i+from] = locDom.m_delx_xi[i+n];
                             m_delx_eta[i+from] = locDom.m_delx_eta[i+n];
                             m_delx_zeta[i+from] = locDom.m_delx_zeta[i+n];
                             m_e[i+from] = locDom.m_e[i+n];
                             m_p[i+from] = locDom.m_p[i+n];
                             m_q[i+from] = locDom.m_q[i+n];
                             m_ql[i+from] = locDom.m_ql[i+n];
                             m_qq[i+from] = locDom.m_qq[i+n];
                             m_v[i+from] = locDom.m_v[i+n];
                             m_volo[i+from] = locDom.m_volo[i+n];
                             m_delv[i+from] = locDom.m_delv[i+n];
                             m_vdov[i+from] = locDom.m_vdov[i+n];
                             m_arealg[i+from] = locDom.m_arealg[i+n];
                             m_ss[i+from] = locDom.m_ss[i+n];
                             m_elemMass[i+from] = locDom.m_elemMass[i+n];                             
                        }
                    }
                }
            }
        }
    }
}


void mapGlobaltoLocalNodes(int numRanks,
                              int myRank,
                              Domain& locDom,
                              std::vector<double>& m_fx,
                              std::vector<double>& m_fy,
                              std::vector<double>& m_fz,
                              std::vector<double>& m_x,
                              std::vector<double>& m_y,
                              std::vector<double>& m_z,
                              std::vector<double>& m_xd,
                              std::vector<double>& m_yd,
                              std::vector<double>& m_zd,
                              std::vector<double>& m_xdd,
                              std::vector<double>& m_ydd,
                              std::vector<double>& m_zdd,
                              std::vector<double>& m_nodalMass)
{
   int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    //Laik_Space* space = laik_partitioning_get_space(p);
    //const Laik_Slice* slice = laik_space_asslice(space);
    //int edgeElems= (int) (cbrt( (domain.sizeX()+1) / numRanks ) + 0.1 );

    int Nx=locDom.sizeX()+1,Ny=locDom.sizeX()+1,Nz=locDom.sizeX()+1;

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


    //printf("Rank: %d, Rx, Ry, Rz = %d, Lz = %d, Ly= %d, Pxy= %d\n", myRank, Rx, Lx, Ly, Pxy);

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
                           locDom.m_fx[i+n] = m_fx[i+from];
                           locDom.m_fy[i+n] = m_fy[i+from];
                           locDom.m_fz[i+n] = m_fz[i+from];
                           locDom.m_nodalMass[i+n] = m_nodalMass[i+from];
                           locDom.m_x[i+n] = m_x[i+from];
                           locDom.m_y[i+n] = m_y[i+from];
                           locDom.m_z[i+n] = m_z[i+from];
                           locDom.m_xd[i+n] = m_xd[i+from];
                           locDom.m_yd[i+n] = m_yd[i+from];
                           locDom.m_zd[i+n] = m_zd[i+from];
                           locDom.m_xdd[i+n] = m_xdd[i+from];
                           locDom.m_xdd[i+n] = m_ydd[i+from];
                           locDom.m_zdd[i+n] = m_zdd[i+from];              
                        }
                    }
                }
            }
        }
    }
}

void mapGlobaltoLocalElements(int numRanks,
                              int myRank,
                              Domain& locDom,
                              std::vector<double>& m_delv_xi,
                              std::vector<double>& m_delv_eta,
                              std::vector<double>& m_delv_zeta,
                              std::vector<double>& m_dxx,
                              std::vector<double>& m_dyy,
                              std::vector<double>& m_dzz,
                              std::vector<double>& m_delx_xi,
                              std::vector<double>& m_delx_eta,
                              std::vector<double>& m_delx_zeta,
                              std::vector<double>& m_e,
                              std::vector<double>& m_p,
                              std::vector<double>& m_q,
                              std::vector<double>& m_ql,
                              std::vector<double>& m_qq,
                              std::vector<double>& m_v,
                              std::vector<double>& m_volo,
                              std::vector<double>& m_delv,
                              std::vector<double>& m_vdov,
                              std::vector<double>& m_arealg,
                              std::vector<double>& m_ss,
                              std::vector<double>& m_elemMass
                             )
{
  
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    //Laik_Space* space = laik_partitioning_get_space(p);
    //const Laik_Slice* slice = laik_space_asslice(space);
    //int edgeElems= (int) (cbrt( (domain.sizeX()+1) / numRanks ) + 0.1 );

    int Nx=locDom.sizeX(),Ny=locDom.sizeX(),Nz=locDom.sizeX();
    int Nxh=locDom.sizeX()+1, Nyh=locDom.sizeX()+1, Nzh=locDom.sizeX()+1;

    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*Nx;
    int Ly = Ry*Ny;
    //int Lz = Rz*Nz;
    int Pxy= Lx*Ly;
    //int Pxz= Lx*Lz;
    //int Pyz= Ly*Lz;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this
    // task
    unsigned long  from, to,n;
    int r=0;
    int nx=0;
    int tag=0;
    for (int rz = 0; rz < Rz; rz++)
    {
        for (int ry = 0; ry < Ry; ry++)
        {
            for (int rx = 0; rx < Rx; rx++)
            {
                r = rx + ry*Rx + rz*Rx*Ry; // task number
                if(r!=myRank) continue;
                // loop over z and x  to create the slices in the
                // partitioning
                for (int ny = 0; ny < Ny; ny++)
                {
                    for (int nz = 0; nz < Nz; nz++)
                    {   
                        // tag = global index where nx = 0 + safety shift = Ny+10
                        nx=0;
                        tag = nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz + Ny*10;
                        nx=0;
                        from = nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx=Nx;
                        to = nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx=0;
                        n = nx + ny*Nx + nz*Nx*Ny;

                        #pragma ivdep
                        for(int i=0; i<(to-from); i++){
                           locDom.m_delv_xi[i+n] =  m_delv_xi[i+from];
                           locDom.m_delv_eta[i+n] =  m_delv_eta[i+from];
                           locDom.m_delv_zeta[i+n] = m_delv_zeta[i+from];
                           locDom.m_dxx[i+n] =  m_dxx[i+from];
                           locDom.m_dyy[i+n] =  m_dyy[i+from];
                           locDom.m_dzz[i+n] =  m_dzz[i+from];
                           locDom.m_delx_xi[i+n] =  m_delx_xi[i+from];
                           locDom.m_delx_eta[i+n] =  m_delx_eta[i+from];
                           locDom.m_delx_zeta[i+n] =  m_delx_zeta[i+from];
                           locDom.m_e[i+n] =  m_e[i+from];
                           locDom.m_p[i+n] =  m_p[i+from];
                           locDom.m_q[i+n] =  m_q[i+from];
                           locDom.m_ql[i+n] =  m_ql[i+from];
                           locDom.m_qq[i+n] =  m_qq[i+from];
                           locDom.m_v[i+n] =  m_v[i+from];
                           locDom.m_volo[i+n] =  m_volo[i+from];
                           locDom.m_delv[i+n] =  m_delv[i+from];
                           locDom.m_vdov[i+n] =  m_vdov[i+from];
                           locDom.m_arealg[i+n] =  m_arealg[i+from];
                           locDom.m_ss[i+n] =  m_ss[i+from];
                           locDom.m_elemMass[i+n] =  m_elemMass[i+from];
                        }
                    }
                }
            }
        }
    }
}




static void PrintCommandLineOptions(char *execname, int myRank)
{
   if (myRank == 0) {

      printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf("\n\n");
   }
}

static void ParseError(const char *message, int myRank)
{
   if (myRank == 0) {
      printf("%s\n", message);
#if USE_MPI      
        assert(-1);
      MPI_Abort(mpi.current_comm, -1);
#else
      exit(-1);
#endif
   }
}

void ParseCommandLineOptions(int argc, char *argv[],
                             Int_t myRank, struct cmdLineOpts *opts)
{
   if(argc > 1) {
      int i = 1;

      while(i < argc) {
         int ok;
         /* -i <iterations> */
         if(strcmp(argv[i], "-i") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -i", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->its));
            if(!ok) {
               ParseError("Parse Error on option -i integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -s <size, sidelength> */
         else if(strcmp(argv[i], "-s") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -s\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->nx));
            if(!ok) {
               ParseError("Parse Error on option -s integer value required after argument\n", myRank);
            }
            i+=2;
         }
	 /* -r <numregions> */
         else if (strcmp(argv[i], "-r") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -r\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->numReg));
            if (!ok) {
               ParseError("Parse Error on option -r integer value required after argument\n", myRank);
            }
            i+=2;
         }
	 /* -f <numfilepieces> */
         else if (strcmp(argv[i], "-f") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -f\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->numFiles));
            if (!ok) {
               ParseError("Parse Error on option -f integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -p */
         else if (strcmp(argv[i], "-p") == 0) {
            opts->showProg = 1;
            i++;
         }
         /* -q */
         else if (strcmp(argv[i], "-q") == 0) {
            opts->quiet = 1;
            i++;
         }
         else if (strcmp(argv[i], "-b") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -b\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->balance));
            if (!ok) {
               ParseError("Parse Error on option -b integer value required after argument\n", myRank);
            }
            i+=2;
         }
         else if (strcmp(argv[i], "-c") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -c\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->cost));
            if (!ok) {
               ParseError("Parse Error on option -c integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -v */
         else if (strcmp(argv[i], "-v") == 0) {
#if VIZ_MESH            
            opts->viz = 1;
#else
            ParseError("Use of -v requires compiling with -DVIZ_MESH\n", myRank);
#endif
            i++;
         }
         else if (strcmp(argv[i], "-ps") == 0){
            i+=2;
         }
         /* -h */
         else if (strcmp(argv[i], "-h") == 0) {
            PrintCommandLineOptions(argv[0], myRank);
#if USE_MPI            
        assert(-1);
            MPI_Abort(mpi.current_comm, 0);
#else
            exit(0);
#endif
         }
         
         else {
            char msg[80];
            PrintCommandLineOptions(argv[0], myRank);
            sprintf(msg, "ERROR: Unknown command line argument: %s\n", argv[i]);
            ParseError(msg, myRank);
         }
      }
   }
}

/////////////////////////////////////////////////////////////////////

void VerifyAndWriteFinalOutput(Real_t elapsed_time,
                               Domain& locDom,
                               Int_t nx,
                               Int_t numRanks)
{
   // GrindTime1 only takes a single domain into account, and is thus a good way to measure
   // processor speed indepdendent of MPI parallelism.
   // GrindTime2 takes into account speedups from MPI parallelism.
   // Cast to 64-bit integer to avoid overflows.
   Int8_t nx8 = nx;
   Real_t grindTime1 = ((elapsed_time*1e6)/locDom.cycle())/(nx8*nx8*nx8);
   Real_t grindTime2 = ((elapsed_time*1e6)/locDom.cycle())/(nx8*nx8*nx8*numRanks);

   Index_t ElemId = 0;
   std::cout << "Run completed:\n";
   std::cout << "   Problem size        =  " << nx       << "\n";
   std::cout << "   MPI tasks           =  " << numRanks << "\n";
   std::cout << "   Iteration count     =  " << locDom.cycle() << "\n";
   std::cout << "   Final Origin Energy =  ";
   std::cout << std::scientific << std::setprecision(6);
   std::cout << std::setw(12) << locDom.e(ElemId) << "\n";

   Real_t   MaxAbsDiff = Real_t(0.0);
   Real_t TotalAbsDiff = Real_t(0.0);
   Real_t   MaxRelDiff = Real_t(0.0);

   for (Index_t j=0; j<nx; ++j) {
      for (Index_t k=j+1; k<nx; ++k) {
         Real_t AbsDiff = FABS(locDom.e(j*nx+k)-locDom.e(k*nx+j));
         TotalAbsDiff  += AbsDiff;

         if (MaxAbsDiff <AbsDiff) MaxAbsDiff = AbsDiff;

         Real_t RelDiff = AbsDiff / locDom.e(k*nx+j);

         if (MaxRelDiff <RelDiff)  MaxRelDiff = RelDiff;
      }
   }

   // Quick symmetry check
   std::cout << "   Testing Plane 0 of Energy Array on rank 0:\n";
   std::cout << "        MaxAbsDiff   = " << std::setw(12) << MaxAbsDiff   << "\n";
   std::cout << "        TotalAbsDiff = " << std::setw(12) << TotalAbsDiff << "\n";
   std::cout << "        MaxRelDiff   = " << std::setw(12) << MaxRelDiff   << "\n";

   // Timing information
   std::cout.unsetf(std::ios_base::floatfield);
   std::cout << std::setprecision(2);
   std::cout << "\nElapsed time         = " << std::setw(10) << elapsed_time << " (s)\n";
   std::cout << std::setprecision(8);
   std::cout << "Grind time (us/z/c)  = "  << std::setw(10) << grindTime1 << " (per dom)  ("
             << std::setw(10) << elapsed_time << " overall)\n";
   std::cout << "FOM                  = " << std::setw(10) << 1000.0/grindTime2 << " (z/s)\n\n";

   return ;
}
