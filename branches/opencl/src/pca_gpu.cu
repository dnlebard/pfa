#include <stdio.h>
#include <stdbool.h>
#include <cuda.h>

/*********************************
HEADER SECTION
//////////////////////////////////
*********************************/
//! Pointers for storing the protein data
float4* prot_dev;  // protein atoms on the device
float4* prot_host; // protein atoms on the host

float4* prot_crys_dev;  // protein crystal structure on the device
float4* prot_crys_host; // protein crystal structure on the host

float4* prot_avg_dev;  // average protein structure on the device
float4* prot_avg_host; // average protein structure on the host

float* prot_vec_dev;    // protein configuration on the device (in vector format)
float* prot_vec_host;   // protein configuration on the host (in vector format)

float* prot_avg_vec_dev;  // protein avg. configuration on the device (in vector format)
float* prot_avg_vec_host; // protein avg. configuration on the host (in vector format)

float* covar_dev;    // covariance matrix on the device
                     // the hosts version will reside in fortran's memory space...good luck!

float* kmat_11_dev;   // the 9 invidivual elements of the 
float* kmat_11_host;  // Kabsch matrix for each protein atom
float* kmat_12_dev;   // on the device and host
float* kmat_12_host;  //
float* kmat_13_dev;   //
float* kmat_13_host;  //

float* kmat_21_dev;   // 
float* kmat_21_host;  // 
float* kmat_22_dev;   // 
float* kmat_22_host;  //
float* kmat_23_dev;   //
float* kmat_23_host;  //

float* kmat_31_dev;   // 
float* kmat_31_host;  // 
float* kmat_32_dev;   // 
float* kmat_32_host;  //
float* kmat_33_dev;   //
float* kmat_33_host;  //
////////////////////////

// System-wide numbers
int NProtein; // the number of protein atoms
int NSelection; // the number of protein atoms
int NProteinDOF; // the number of degrees of freedom for all protein atoms
int NSelectionDOF; // the number of degrees of freedom for all selected atoms
int iPStart; // the starting index of the protein (0-based)
int iPEnd; // the starting index of the protein (0-based)

// Size infos
unsigned int memProteinF4;
unsigned int memProteinDOF;
unsigned int memKabschF1;
unsigned int memSelectionF4;
unsigned int memSelectionDOF;
unsigned int memKabschSelF1;
//////////////////////////////////
// Simple utility function to check for CUDA runtime errors
void checkCUDAError22(const char* msg);

extern "C" { extern "C" void alloc_gpu_pca_(int* rank, int* iPrStart, int* iPrEnd, float* qrms, bool* doSelect, int* selIndexes, int* iSelect, int* iGPUs, bool* doGPUInit);}

extern "C" { extern "C" void dealloc_gpu_pca_();}

extern "C" {extern "C" void get_avg_gpu_(float *xavg,float *yavg,float *zavg);}

extern "C" {extern "C" void get_avg_sel_gpu_(float *xavg,float *yavg,float *zavg);}

extern "C" {extern "C" void init_covar_gpu_();}

extern "C" {extern "C" void init_covar_sel_gpu_();}

extern "C" {extern "C" void accum_covar_gpu_(float *xyz, float* qrms, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3);}

extern "C" {extern "C" void accum_avg_gpu_(float *xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3);}

extern "C" {extern "C" void accum_avg_sel_gpu_(float *xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3);}

extern "C" {extern "C" void store_crys_gpu_(float *xyz, float* qrms,bool* setVectors);}

extern "C" {extern "C" void store_crys_sel_gpu_(float *xyz, float* qrms,bool* setVectors,int* selIndexes);}

extern "C" {extern "C" void get_crys_gpu_(float *xyz);}

extern "C" {extern "C" void get_covar_gpu_(float *covar_cpu);}

extern "C" {extern "C" void get_covar_sel_gpu_(float *covar_cpu);}

extern "C" {extern "C" void calc_kmat_gpu_(float* xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3);}

extern "C" {extern "C" void calc_kmat_sel_gpu_(float* xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3);}

extern "C" {extern "C" void accum_covar_sel_gpu_(float *xyz, float* qrms, int* selIndexes, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3);}
/*********************************
END HEADER SECTION
//////////////////////////////////
*********************************/
/*
    initCovar_GPU<<<gridSize, blockSize>>>(covar_dev);
*/
__global__ void initCovar_GPU(int nDOF,float* covar_dev)
{

     // 
     unsigned int idy_covar = blockIdx.x*blockDim.x + threadIdx.x;

     if (idy_covar < nDOF)
     {// then we actually do something...

       // loop over one column of the covariance matrix
       for(int idx_covar=0;idx_covar<nDOF;++idx_covar)
       {
        // init to 0
        covar_dev[idy_covar*nDOF + idx_covar] = 0.0f;
        // init to 0

       }

     }

   __syncthreads();
   return; 
}
/*
    accumCovar_GPU<<<gridSize, blockSize>>>(NProteinDOF,prot_vec_dev,prot_avg_vec_dev,covar_dev);
*/
__global__ void accumCovar_GPU(int nDOF,float* prot_vec_dev,float* prot_avg_vec_dev,float* covar_dev)
{

     // 
     unsigned int idy_covar = blockIdx.x*blockDim.x + threadIdx.x;

     if (idy_covar < nDOF)
     {// then we actually do something...

       float avg_y = prot_avg_vec_dev[idy_covar];
       float vec_y = prot_vec_dev[idy_covar];
       float dy    = vec_y - avg_y;

       // loop over one column of the covariance matrix
       for(int idx_covar=0;idx_covar<nDOF;++idx_covar)
       {
        float avg_x = prot_avg_vec_dev[idx_covar];
        float vec_x = prot_vec_dev[idx_covar];
        float dx    = vec_x - avg_x;

        float oldVal = covar_dev[idy_covar*nDOF + idx_covar];

        float cVal  = dx*dy;
              cVal += oldVal;
//        covar_dev[idy_covar*nDOF + idx_covar] += cVal;

//         covar_dev[idy_covar*nDOF + idx_covar] += (dx*dy);

        // storing data in column major order like fortran
          covar_dev[idy_covar*nDOF + idx_covar] = cVal;
//        covar_dev[idx_covar*nDOF + idy_covar] = cVal;

//        if((idy_covar*nDOF + idx_covar) == 13638248)
//        {
//          printf("accumCovar(DEV) covar_dev(%i,%i) = %f \n",idx_covar,idy_covar,cVal);
//          printf("                nDOF = %i, nDOF^2 = %i \n",nDOF,nDOF*nDOF);
//          printf("                cVal = %f, oldVal = %f \n",dx*dy,oldVal);
//          printf("                vec(i) = %f, avg(i) = %f \n",vec_x,avg_x);
//          printf("                vec(j) = %f, avg(j) = %f \n",vec_y,avg_y);
//         }

       }


     }

   __syncthreads();
   return; 
}


/*
     
rotateProt_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtein,prot_dev,rmat11, rmat12, rmat13, rmat21, rmat22, rmat23, rmat31, rmat32, rmat33);

*/
__global__ void rotateProt_GPU(int nP,float4* prot_dev, float rmat11, float rmat12, float rmat13, float rmat21, float rmat22, float rmat23, float rmat31, float rmat32, float rmat33)
{
     unsigned int idx_prot = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_prot < nP)
     {// then we actually do something...

  
      // get our protein atom's position
      float4 prot_pos = prot_dev[idx_prot];//tex1Dfetch(water_tex, idx_water);   
      float pX = prot_pos.x;
      float pY = prot_pos.y;
      float pZ = prot_pos.z;

      // rotate the protein atom
      float x  = rmat11*pX;
      x += rmat12*pY;
      x += rmat13*pZ;

      float y  = rmat21*pX;
      y += rmat22*pY;
      y += rmat23*pZ;

      float z  = rmat31*pX;
      z += rmat32*pY;
      z += rmat33*pZ;


       float4 prot_pos_new;
       prot_pos_new.x = x;
       prot_pos_new.y = y;
       prot_pos_new.z = z;

        prot_dev[idx_prot] = prot_pos_new;

//      prot_crys_dev[idx_prot].x = pX;
//      prot_crys_dev[idx_prot].y = pY;
//      prot_crys_dev[idx_prot].z = pZ;

/*
      if(fabs(pX) > 1000000.0 || fabs(pY) > 1000000.0 || fabs(pZ) > 1000000.0)
      {
          printf("rotateProt(DEV) : px = %f, py = %f, pz = %f\n",pX,pY,pZ);
          printf("                r11 = %f r12 = %f r13 = %f \n",rmat11,rmat12,rmat13);
          printf("                r21 = %f r22 = %f r23 = %f \n",rmat21,rmat22,rmat23);
          printf("                r31 = %f r32 = %f r33 = %f \n",rmat31,rmat32,rmat33);
          printf("                prx = %f, pry = %f, prz = %f\n",prot_pos.x,prot_pos.y,prot_pos.z);
      }
*/

     }              
   __syncthreads();
     return; 

}

/*
     
  accumProtStruc_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtein,prot_dev,prot_avg_dev, rmat11, rmat12, rmat13, rmat21, rmat22, rmat23, rmat31, rmat32, rmat33)

*/
__global__ void accumProtStruc_GPU(int nP,float4* prot_dev,float4* prot_avg_dev,float rmat11, float rmat12, float rmat13, float rmat21, float rmat22, float rmat23, float rmat31, float rmat32, float rmat33)
{
     unsigned int idx_prot = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_prot < nP)
     {// then we actually do something...

  
      // get our protein atom's position
      float4 prot_pos = prot_dev[idx_prot];//tex1Dfetch(water_tex, idx_water);   
      float pX = prot_pos.x;
      float pY = prot_pos.y;
      float pZ = prot_pos.z;

      // rotate the protein atom
      float x  = rmat11*pX;
      x += rmat12*pY;
      x += rmat13*pZ;

      float y  = rmat21*pX;
      y += rmat22*pY;
      y += rmat23*pZ;

      float z  = rmat31*pX;
      z += rmat32*pY;
      z += rmat33*pZ;

      prot_avg_dev[idx_prot].x += x;
      prot_avg_dev[idx_prot].y += y;
      prot_avg_dev[idx_prot].z += z;


     }              
   __syncthreads();
     return; 

}

/*
     
   calckMat_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtein,prot_dev,prot_avg_dev,kmat_11_dev, kmat_12_dev, kmat_13_dev, kmat_21_dev, kmat_22_dev, kmat_23_dev, kmat_31_dev, kmat_32_dev, kmat_33_dev)

*/
__global__ void calckMat_GPU(int nP,float4* prot_dev,float4* prot_crys_dev, float* kmat_11_dev, float* kmat_12_dev, float* kmat_13_dev, float* kmat_21_dev, float* kmat_22_dev, float* kmat_23_dev, float* kmat_31_dev, float* kmat_32_dev, float* kmat_33_dev)
{
     unsigned int idx_prot = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_prot < nP)
     {// then we actually do something...

  
      // get our protein atom's position
      float4 prot_pos = prot_dev[idx_prot];//tex1Dfetch(water_tex, idx_water);   
      float pX = prot_pos.x;
      float pY = prot_pos.y;
      float pZ = prot_pos.z;

      // get our crystal atom's position
      float4 prot_crys = prot_crys_dev[idx_prot];//tex1Dfetch(water_tex, idx_water);   
      float pcX = prot_crys.x;
      float pcY = prot_crys.y;
      float pcZ = prot_crys.z;
      float mass = prot_crys.w; // w position has the mass

      float k11;
      float k12;
      float k13;

      float k21;
      float k22;
      float k23;

      float k31;
      float k32;
      float k33;

      k11 = mass * pcX * pX;
      k12 = mass * pcX * pY;
      k13 = mass * pcX * pZ;

      k21 = mass * pcY * pX;
      k22 = mass * pcY * pY;
      k23 = mass * pcY * pZ;

      k31 = mass * pcZ * pX;
      k32 = mass * pcZ * pY;
      k33 = mass * pcZ * pZ;

      kmat_11_dev[idx_prot] = k11;
      kmat_12_dev[idx_prot] = k12;
      kmat_13_dev[idx_prot] = k13;

      kmat_21_dev[idx_prot] = k21;
      kmat_22_dev[idx_prot] = k22;
      kmat_23_dev[idx_prot] = k23;

      kmat_31_dev[idx_prot] = k31;
      kmat_32_dev[idx_prot] = k32;
      kmat_33_dev[idx_prot] = k33;

     }              
   __syncthreads();
     return; 

}
/////////////////////////////////////////////////////
////////////// END GPU KERKEL SECTION ///////////////
/////////////////////////////////////////////////////
extern "C" void alloc_gpu_pca_(int* rank, int* iPrStart, int* iPrEnd, float* qrms, bool* doSelect, int* selIndexes, int* iSelect, int* iGPUs, bool* doGPUInit)
{

   iPStart = *iPrStart - 1;
   NProtein = (*iPrEnd - *iPrStart) + 1;
   NProteinDOF = NProtein * 3;
   iPEnd = NProteinDOF - 3; // the ending protein atomic index
   NSelection = *iSelect;
   NSelectionDOF = NSelection * 3;

   //char *methodName;
   //methodName = "alloc_gpu_pca";
   //alloc_gpu_rank_(rank, iGPUs, methodName);

   if(*doGPUInit)
   {
   // Find our local GPU rank on the node
    char cudaOut[100];
    int iChr;
    int rank_gpu = *rank % *iGPUs;
    cudaSetDevice(rank_gpu); // THIS COULD BE A DANGEROUS PLACE TO MAKE THIS CALL
    iChr = sprintf(cudaOut,"cudaSetDevice -- setting device to %i on rank %i \n",rank_gpu,*rank);
    if(iChr < 0 ){printf("SHOULD I CARE?\n");}
    checkCUDAError22(cudaOut);
    *doGPUInit = false;
   }
   //cudaSetDevice(*rank); // THIS COULD BE A DANGEROUS PLACE TO MAKE THIS CALL


   //------------------------------------------
   // PUT THIS IN AN ALLOCATION ROUTINE
   //------------------------------------------
   memSelectionF4 = sizeof(float4) *  NSelection;
   memKabschSelF1  = sizeof(float) * NSelection;
   memSelectionDOF = sizeof(float) * NSelectionDOF;

   memProteinF4 = sizeof(float4) *  NProtein;
   memKabschF1  = sizeof(float) * NProtein;
   memProteinDOF = sizeof(float) * NProteinDOF;
   if(*doSelect)
   {

     kmat_11_host = (float*)malloc(memKabschSelF1);
     kmat_12_host = (float*)malloc(memKabschSelF1);
     kmat_13_host = (float*)malloc(memKabschSelF1);

     kmat_21_host = (float*)malloc(memKabschSelF1);
     kmat_22_host = (float*)malloc(memKabschSelF1);
     kmat_23_host = (float*)malloc(memKabschSelF1);

     kmat_31_host = (float*)malloc(memKabschSelF1);
     kmat_32_host = (float*)malloc(memKabschSelF1);
     kmat_33_host = (float*)malloc(memKabschSelF1);

     prot_host   = (float4*)malloc(memSelectionF4);
     prot_avg_host   = (float4*)malloc(memSelectionF4);
     prot_crys_host   = (float4*)malloc(memSelectionF4);
     prot_vec_host = (float*)malloc(memSelectionDOF);
     prot_avg_vec_host = (float*)malloc(memSelectionDOF);

     cudaMalloc( (void **) &kmat_11_dev,memKabschSelF1);
     cudaMalloc( (void **) &kmat_12_dev,memKabschSelF1);
     cudaMalloc( (void **) &kmat_13_dev,memKabschSelF1);


     cudaMalloc( (void **) &kmat_21_dev,memKabschSelF1);
     cudaMalloc( (void **) &kmat_22_dev,memKabschSelF1);
     cudaMalloc( (void **) &kmat_23_dev,memKabschSelF1);

     cudaMalloc( (void **) &kmat_31_dev,memKabschSelF1);
     cudaMalloc( (void **) &kmat_32_dev,memKabschSelF1);
     cudaMalloc( (void **) &kmat_33_dev,memKabschSelF1);

     cudaMalloc( (void **) &prot_dev,memSelectionF4);
     cudaMalloc( (void **) &prot_avg_dev,memSelectionF4);
     cudaMalloc( (void **) &prot_crys_dev,memSelectionF4);
     cudaMalloc( (void **) &prot_vec_dev,memSelectionDOF);
     cudaMalloc( (void **) &prot_avg_vec_dev,memSelectionDOF);

     cudaMalloc( (void **) &covar_dev,memSelectionDOF * NSelectionDOF);

     // now, copy over all the protein masses to the gpu
     float* protMass;
     protMass = &qrms[2]; /// points to 2 for the mass position

     int j,k;
     // loop over all protein atoms and fill up the protein array
     for(int i = 0;i < NSelection;i++) 
     {
      k=selIndexes[i];

      j = 3*k - 3;
      // fill up the proteins masses on the host in the w position of the float4s
      prot_avg_host[i].x = 0.0f;
      prot_avg_host[i].y = 0.0f;
      prot_avg_host[i].z = 0.0f;
      prot_avg_host[i].w = protMass[j];
     }
      // copy the masses to the device via the average structure
     cudaMemcpy(prot_avg_dev, prot_avg_host, memSelectionF4, cudaMemcpyHostToDevice);

   } else {
     // the full protein (non-selection)
     kmat_11_host = (float*)malloc(memKabschF1);
     kmat_12_host = (float*)malloc(memKabschF1);
     kmat_13_host = (float*)malloc(memKabschF1);

     kmat_21_host = (float*)malloc(memKabschF1);
     kmat_22_host = (float*)malloc(memKabschF1);
     kmat_23_host = (float*)malloc(memKabschF1);

     kmat_31_host = (float*)malloc(memKabschF1);
     kmat_32_host = (float*)malloc(memKabschF1);
     kmat_33_host = (float*)malloc(memKabschF1);

     prot_host   = (float4*)malloc(memProteinF4);
     prot_avg_host   = (float4*)malloc(memProteinF4);
     prot_crys_host   = (float4*)malloc(memProteinF4);
     prot_vec_host = (float*)malloc(memProteinDOF);
     prot_avg_vec_host = (float*)malloc(memProteinDOF);

     cudaMalloc( (void **) &kmat_11_dev,memKabschF1);
     cudaMalloc( (void **) &kmat_12_dev,memKabschF1);
     cudaMalloc( (void **) &kmat_13_dev,memKabschF1);


     cudaMalloc( (void **) &kmat_21_dev,memKabschF1);
     cudaMalloc( (void **) &kmat_22_dev,memKabschF1);
     cudaMalloc( (void **) &kmat_23_dev,memKabschF1);

     cudaMalloc( (void **) &kmat_31_dev,memKabschF1);
     cudaMalloc( (void **) &kmat_32_dev,memKabschF1);
     cudaMalloc( (void **) &kmat_33_dev,memKabschF1);

     cudaMalloc( (void **) &prot_dev,memProteinF4);
     cudaMalloc( (void **) &prot_avg_dev,memProteinF4);
     cudaMalloc( (void **) &prot_crys_dev,memProteinF4);
     cudaMalloc( (void **) &prot_vec_dev,memProteinDOF);
     cudaMalloc( (void **) &prot_avg_vec_dev,memProteinDOF);

     cudaMalloc( (void **) &covar_dev,memProteinDOF * NProteinDOF);

     // now, copy over all the protein masses to the gpu
     float* protMass;
     protMass = &qrms[2]; /// points to 2 for the mass position

     int j = -3;
     // loop over all protein atoms and fill up the protein array
     for(int i = 0;i < NProtein;i++) 
     {
      j+=3;
      // fill up the proteins masses on the host in the w position of the float4s
      prot_avg_host[i].x = 0.0f;
      prot_avg_host[i].y = 0.0f;
      prot_avg_host[i].z = 0.0f;
      prot_avg_host[i].w = protMass[j];
     }
      // copy the masses to the device via the average structure
     cudaMemcpy(prot_avg_dev, prot_avg_host, memProteinF4, cudaMemcpyHostToDevice);
  
   } // end test for selection
}

extern "C" void dealloc_gpu_pca_()
{
   //------------------------------------------
   // PUT THIS IN A DEALLOCATION ROUTINE
   //------------------------------------------
   // free device memory

    cudaFree(kmat_11_dev);
    cudaFree(kmat_12_dev);
    cudaFree(kmat_13_dev);

    cudaFree(kmat_21_dev);
    cudaFree(kmat_22_dev);
    cudaFree(kmat_23_dev);

    cudaFree(kmat_31_dev);
    cudaFree(kmat_32_dev);
    cudaFree(kmat_33_dev);

    cudaFree(prot_dev);
    cudaFree(prot_avg_dev);
    cudaFree(prot_crys_dev);
    cudaFree(prot_vec_dev);
    cudaFree(prot_avg_vec_dev);
    cudaFree(covar_dev);

   // free up host memory     

      delete[] kmat_11_host;
      kmat_11_host = NULL;
      delete[] kmat_12_host;
      kmat_12_host = NULL;
      delete[] kmat_13_host;
      kmat_13_host = NULL;

      delete[] kmat_21_host;
      kmat_21_host = NULL;
      delete[] kmat_22_host;
      kmat_22_host = NULL;
      delete[] kmat_23_host;
      kmat_23_host = NULL;

      delete[] kmat_31_host;
      kmat_31_host = NULL;
      delete[] kmat_32_host;
      kmat_32_host = NULL;
      delete[] kmat_33_host;
      kmat_33_host = NULL;

      delete[] prot_host;
      prot_host = NULL;

      delete[] prot_avg_host;
      prot_avg_host = NULL;

      delete[] prot_crys_host;
      prot_crys_host = NULL;

      delete[] prot_vec_host;
      prot_crys_host = NULL;

      delete[] prot_avg_vec_host;
      prot_avg_vec_host = NULL;
   //------------------------------------------
   // END DEALLOCATION ROUTINE
   //------------------------------------------

}

/*

init_covar_gpu_: (called from fortran as "init_covar_gpu") is the fortran interface function for the 
                kernel to init the proteins covariance matrix to 0.0f

*/
extern "C" void init_covar_gpu_()
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = NProteinDOF / numThreadsPerBlock + (NProteinDOF % numThreadsPerBlock == 0?0:1);

//   printf("NOW INITING THE COVARIANCE MATRIX ON GPU...\n");
   // now accumulate the covariance matrix
   initCovar_GPU<<<numBlocks, numThreadsPerBlock>>>(NProteinDOF,covar_dev);

}

/*

init_covar_sel_gpu_: (called from fortran as "init_covar_sel_gpu") is the fortran interface function for the 
                kernel to init the selection atoms covariance matrix to 0.0f

*/
extern "C" void init_covar_sel_gpu_()
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock;
   if(NSelectionDOF >= 1000) {
     numThreadsPerBlock = 32;
   } else {
     numThreadsPerBlock = 256;
   }

   int numBlocks = NSelectionDOF / numThreadsPerBlock + (NSelectionDOF % numThreadsPerBlock == 0?0:1);

//   printf("NOW INITING THE COVARIANCE MATRIX ON GPU...\n");
   // now accumulate the covariance matrix
   initCovar_GPU<<<numBlocks, numThreadsPerBlock>>>(NSelectionDOF,covar_dev);

}


/*

accum_covar_gpu_: (called from fortran as "accum_covar_gpu") is the fortran interface function for the 
                kernel to accumulate the proteins covariance matrix on the GPU after rotating the 
                protein to fit the crystal structure and translating it back in the central box

*/
extern "C" void accum_covar_gpu_(float *xyz, float* qrms, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocksDOF = NProteinDOF / numThreadsPerBlock + (NProteinDOF % numThreadsPerBlock == 0?0:1);
   int numBlocks = NProtein / numThreadsPerBlock + (NProtein % numThreadsPerBlock == 0?0:1);

//   printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
//   printf("RUNNING ON GPU WITH %i NProteinDOF\n",NProteinDOF);
   
   float* xx;
   float* yy;
   float* zz;
   float rmat11;
   float rmat12;
   float rmat13;
   float rmat21;
   float rmat22;
   float rmat23;
   float rmat31;
   float rmat32;
   float rmat33;

   double r11 = *one1;
   double r12 = *one2;
   double r13 = *one3;
   double r21 = *two1;
   double r22 = *two2;
   double r23 = *two3;
   double r31 = *three1;
   double r32 = *three2;
   double r33 = *three3;

   rmat11 = (float)r11;
   rmat12 = (float)r12;
   rmat13 = (float)r13;
   rmat21 = (float)r21;
   rmat22 = (float)r22;
   rmat23 = (float)r23;
   rmat31 = (float)r31;
   rmat32 = (float)r32;
   rmat33 = (float)r33;

//   printf("IN accumCovar(HOST) --> 1st CALL\n");
//   printf("Rmat (%f %f %f)\n",rmat11,rmat12,rmat13);
//   printf("Rmat (%f %f %f)\n",rmat21,rmat22,rmat23);
//   printf("Rmat (%f %f %f)\n",rmat31,rmat32,rmat33);
/*
   rmat11 = (float)*one1;
   rmat12 = (float)*one2;
   rmat13 = (float)*one3;
   rmat21 = (float)*two1;
   rmat22 = (float)*two2;
   rmat23 = (float)*two3;
   rmat31 = (float)*three1;
   rmat32 = (float)*three2;
   rmat33 = (float)*three3;
*/
/*
   xx = &xyz[iPStart+0];
   yy = &xyz[iPStart+1];  THIS WAS REMOVED BECAUSE nodeConf shouldnt be given as input
   zz = &xyz[iPStart+2];
*/
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];


    int j = -3;
    for(int i = 0;i < NProtein;i++) 
   {
    j+=3;
    // fill up the host with the protein position at the current timestep
    prot_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);

//     if(prot_host[i].x > 100000000 || prot_host[i].y > 100000000 || prot_host[i].z > 100000000)
//     {  
//          printf("accumCovar(HOST) prot_host[i].x = %f, prot_host[i].y = %f,prot_host[i].z = %f \n",prot_host[i].x,prot_host[i].y,prot_host[i].z);
//     }
   }

    // copy the current protein position to the device
   cudaMemcpy(prot_dev, prot_host, memProteinF4, cudaMemcpyHostToDevice);

//   printf("ROTATING THE PROTEIN ON GPU...\n");

   // rotate the protein on the GPU
   rotateProt_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtein,prot_dev,rmat11, rmat12, rmat13, rmat21, rmat22, rmat23, rmat31, rmat32, rmat33);
//   printf("ROTATED\n");

//   printf("COPYING THE ROTATED PROTEIN BACK TO THE GPU...\n");
    // copy the rotated protein back to the host
   cudaMemcpy(prot_host, prot_dev, memProteinF4, cudaMemcpyDeviceToHost);
//   printf("COPYIED.\n");

   j = -3;
   int k = -1;

   float* protMass = &qrms[2]; /// points to 2 for the mass position
   float  sqrtMass;

//   printf("ABOUT TO LOOP, NProtein = %i.\n",NProtein);
    for(int i = 0;i < NProtein;++i) 
   {
    j+=3;

    sqrtMass = sqrtf(protMass[j]);

    // fill up the host with the protein position at the current timestep
    k++;
    prot_vec_host[k] = prot_host[i].x*sqrtMass;
    k++;
    prot_vec_host[k] = prot_host[i].y*sqrtMass;
    k++;
    prot_vec_host[k] = prot_host[i].z*sqrtMass;

//     if(i == (NProtein - 1))
//     {  
//          printf("accumCovar(HOST) k = %i, k-1 = %i,k-2 = %i \n",k, k-1,k-2);
//          printf("                prot_host[j].x,y,z = %f,%f,%f  \n",prot_host[i].x,prot_host[i].y,prot_host[i].z);
//          printf("                sqrtMass = %f \n",sqrtMass);
//     }

   }

    // copy the current protein position to the device
   cudaMemcpy(prot_vec_dev, prot_vec_host, memProteinDOF, cudaMemcpyHostToDevice);

//   printf("NOW ACCUMULATING THE COVARIANCE MATRIX ON GPU...\n");
   // now accumulate the covariance matrix
   accumCovar_GPU<<<numBlocksDOF, numThreadsPerBlock>>>(NProteinDOF,prot_vec_dev,prot_avg_vec_dev,covar_dev);

}

/*

accum_covar_sel_gpu_: (called from fortran as "accum_covar_sel_gpu") is the fortran interface function for the 
                kernel to accumulate the selection atoms covariance matrix on the GPU after rotating the 
                selection to fit the crystal structure and translating it back in the central box

*/
extern "C" void accum_covar_sel_gpu_(float *xyz, float* qrms, int* selIndexes, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock;
   if(NSelection >= 1000) {
     numThreadsPerBlock = 32;
   } else {
     numThreadsPerBlock = 256;
   }

   int numBlocksDOF = NSelectionDOF / numThreadsPerBlock + (NSelectionDOF % numThreadsPerBlock == 0?0:1);
   int numBlocks = NSelection / numThreadsPerBlock + (NSelection % numThreadsPerBlock == 0?0:1);

//   printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
//   printf("RUNNING ON GPU WITH %i NProteinDOF\n",NProteinDOF);
   
   float* xx;
   float* yy;
   float* zz;
   float rmat11;
   float rmat12;
   float rmat13;
   float rmat21;
   float rmat22;
   float rmat23;
   float rmat31;
   float rmat32;
   float rmat33;

   double r11 = *one1;
   double r12 = *one2;
   double r13 = *one3;
   double r21 = *two1;
   double r22 = *two2;
   double r23 = *two3;
   double r31 = *three1;
   double r32 = *three2;
   double r33 = *three3;

   rmat11 = (float)r11;
   rmat12 = (float)r12;
   rmat13 = (float)r13;
   rmat21 = (float)r21;
   rmat22 = (float)r22;
   rmat23 = (float)r23;
   rmat31 = (float)r31;
   rmat32 = (float)r32;
   rmat33 = (float)r33;

//   printf("IN accumCovar(HOST) --> 1st CALL\n");
//   printf("Rmat (%f %f %f)\n",rmat11,rmat12,rmat13);
//   printf("Rmat (%f %f %f)\n",rmat21,rmat22,rmat23);
//   printf("Rmat (%f %f %f)\n",rmat31,rmat32,rmat33);
/*
   rmat11 = (float)*one1;
   rmat12 = (float)*one2;
   rmat13 = (float)*one3;
   rmat21 = (float)*two1;
   rmat22 = (float)*two2;
   rmat23 = (float)*two3;
   rmat31 = (float)*three1;
   rmat32 = (float)*three2;
   rmat33 = (float)*three3;
*/
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];
 
    int j = -3;
    for(int i = 0;i < NSelection;i++) 
   {
    j+=3;
    // fill up the host with the protein position at the current timestep
    prot_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);

//     if(prot_host[i].x > 100000000 || prot_host[i].y > 100000000 || prot_host[i].z > 100000000)
//     {  
//          printf("accumCovar(HOST) prot_host[i].x = %f, prot_host[i].y = %f,prot_host[i].z = %f \n",prot_host[i].x,prot_host[i].y,prot_host[i].z);
//     }
   }

    // copy the current protein position to the device
   cudaMemcpy(prot_dev, prot_host, memSelectionF4, cudaMemcpyHostToDevice);

//   printf("ROTATING THE PROTEIN ON GPU...\n");

   // rotate the protein on the GPU
   rotateProt_GPU<<<numBlocks, numThreadsPerBlock>>>(NSelection,prot_dev,rmat11, rmat12, rmat13, rmat21, rmat22, rmat23, rmat31, rmat32, rmat33);
//   printf("ROTATED\n");

//   printf("COPYING THE ROTATED PROTEIN BACK TO THE GPU...\n");
    // copy the rotated protein back to the host
   cudaMemcpy(prot_host, prot_dev, memSelectionF4, cudaMemcpyDeviceToHost);
//   printf("COPYIED.\n");

   int k = -1;

   float* protMass = &qrms[2]; /// points to 2 for the mass position
   float  sqrtMass;

//   printf("ABOUT TO LOOP, NProtein = %i.\n",NProtein);
    for(int i = 0;i < NSelection;++i) 
   {
    j = selIndexes[i];
    j = 3*j - 3;
    sqrtMass = sqrtf(protMass[j]);

    // fill up the host with the protein position at the current timestep
    k++;
    prot_vec_host[k] = prot_host[i].x*sqrtMass;
    k++;
    prot_vec_host[k] = prot_host[i].y*sqrtMass;
    k++;
    prot_vec_host[k] = prot_host[i].z*sqrtMass;

//     if(i == (NProtein - 1))
//     {  
//          printf("accumCovar(HOST) k = %i, k-1 = %i,k-2 = %i \n",k, k-1,k-2);
//          printf("                prot_host[j].x,y,z = %f,%f,%f  \n",prot_host[i].x,prot_host[i].y,prot_host[i].z);
//          printf("                sqrtMass = %f \n",sqrtMass);
//     }

   }

    // copy the current protein position to the device
   cudaMemcpy(prot_vec_dev, prot_vec_host, memSelectionDOF, cudaMemcpyHostToDevice);

//   printf("NOW ACCUMULATING THE COVARIANCE MATRIX ON GPU...\n");
   // now accumulate the covariance matrix
   accumCovar_GPU<<<numBlocksDOF, numThreadsPerBlock>>>(NSelectionDOF,prot_vec_dev,prot_avg_vec_dev,covar_dev);

}


/*

accum_avg_sel_gpu_: (called from fortran as "accum_avg_sel_gpu") is the fortran interface function for the 
                kernel to accumulate the proteins average structure on the GPU after rotating the 
                protein to fit the crystal structure and translating it back in the central box

*/
extern "C" void accum_avg_sel_gpu_(float *xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock;
   if(NSelection >= 1000) {
     numThreadsPerBlock = 32;
   } else {
     numThreadsPerBlock = 256;
   }

   int numBlocks = NSelection / numThreadsPerBlock + (NSelection % numThreadsPerBlock == 0?0:1);

 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtein);
   
   float* xx;
   float* yy;
   float* zz;
   float rmat11;
   float rmat12;
   float rmat13;
   float rmat21;
   float rmat22;
   float rmat23;
   float rmat31;
   float rmat32;
   float rmat33;

   rmat11 = (float)*one1;
   rmat12 = (float)*one2;
   rmat13 = (float)*one3;
   rmat21 = (float)*two1;
   rmat22 = (float)*two2;
   rmat23 = (float)*two3;
   rmat31 = (float)*three1;
   rmat32 = (float)*three2;
   rmat33 = (float)*three3;

   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];
 
    int j = -3;
    for(int i = 0;i < NSelection;i++) 
   {
    j+=3;
    // fill up the host with the protein position at the current timestep
    prot_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);
   }

    // copy the current protein position to the device
   cudaMemcpy(prot_dev, prot_host, memSelectionF4, cudaMemcpyHostToDevice);

   // run the average structure calculation on the GPU
   accumProtStruc_GPU<<<numBlocks, numThreadsPerBlock>>>(NSelection, prot_dev,prot_avg_dev, rmat11,rmat12,rmat13,rmat21,rmat22,rmat23,rmat31,rmat32,rmat33);

}


/*

accum_avg_gpu_: (called from fortran as "accum_avg_gpu") is the fortran interface function for the 
                kernel to accumulate the proteins average structure on the GPU after rotating the 
                protein to fit the crystal structure and translating it back in the central box

*/
extern "C" void accum_avg_gpu_(float *xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = NProtein / numThreadsPerBlock + (NProtein % numThreadsPerBlock == 0?0:1);

 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtein);
   
   float* xx;
   float* yy;
   float* zz;
   float rmat11;
   float rmat12;
   float rmat13;
   float rmat21;
   float rmat22;
   float rmat23;
   float rmat31;
   float rmat32;
   float rmat33;

   rmat11 = (float)*one1;
   rmat12 = (float)*one2;
   rmat13 = (float)*one3;
   rmat21 = (float)*two1;
   rmat22 = (float)*two2;
   rmat23 = (float)*two3;
   rmat31 = (float)*three1;
   rmat32 = (float)*three2;
   rmat33 = (float)*three3;

/*
   xx = &xyz[iPStart+0];
   yy = &xyz[iPStart+1];
   zz = &xyz[iPStart+2];

   here, i think iPStart should be removed because
   we are operating on the proteinCOM structure and not the
   full nodeConf structure
*/

   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];
 
    int j = -3;
    for(int i = 0;i < NProtein;i++) 
   {
    j+=3;
    // fill up the host with the protein position at the current timestep
    prot_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);
   }

    // copy the current protein position to the device
   cudaMemcpy(prot_dev, prot_host, memProteinF4, cudaMemcpyHostToDevice);

   // run the average structure calculation on the GPU
   accumProtStruc_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtein,prot_dev,prot_avg_dev,rmat11, rmat12, rmat13, rmat21, rmat22, rmat23, rmat31, rmat32, rmat33);


}

/*

calc_kmat_gpu_: (called from fortran as "calc_kmat_gpu") is the fortran interface function for the 
                calculation of the 9 elements of the Kabsch matrix.  This routine assumes the crystal
                structure that we will test against has already been copied to the GPU.

*/
extern "C" void calc_kmat_gpu_(float* xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = NProtein / numThreadsPerBlock + (NProtein % numThreadsPerBlock == 0?0:1);

   float* xx;
   float* yy;
   float* zz;

/*
   xx = &xyz[0+iPStart];
   yy = &xyz[1+iPStart];
   zz = &xyz[2+iPStart];

  here, again iPStart should be 0 because we are working on the proteinCOM structure, 
  not the full nodeConf structure.
*/

   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];

    int j = -3;
    for(int i = 0;i < NProtein;i++) 
   {
    j+=3;
    // fill up the host with the protein position at the current timestep
    prot_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);
   }

    // copy the current protein position to the device
   cudaMemcpy(prot_dev, prot_host, memProteinF4, cudaMemcpyHostToDevice);


   // calculate the individual Kabsch elements on the GPU
   calckMat_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtein,prot_dev,prot_crys_dev,kmat_11_dev, kmat_12_dev, kmat_13_dev, kmat_21_dev, kmat_22_dev, kmat_23_dev, kmat_31_dev, kmat_32_dev, kmat_33_dev);

   // copy the element arrays of the Kabsch matrix to the host
   cudaMemcpy(kmat_11_host, kmat_11_dev, memKabschF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_12_host, kmat_12_dev, memKabschF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_13_host, kmat_13_dev, memKabschF1, cudaMemcpyDeviceToHost);

   cudaMemcpy(kmat_21_host, kmat_21_dev, memKabschF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_22_host, kmat_22_dev, memKabschF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_23_host, kmat_23_dev, memKabschF1, cudaMemcpyDeviceToHost);

   cudaMemcpy(kmat_31_host, kmat_31_dev, memKabschF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_32_host, kmat_32_dev, memKabschF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_33_host, kmat_33_dev, memKabschF1, cudaMemcpyDeviceToHost);

    double k11, k12, k13;
    double k21, k22, k23;
    double k31, k32, k33;

    k11 = 0.00;
    k12 = 0.00;
    k13 = 0.00;

    k21 = 0.00;
    k22 = 0.00;
    k23 = 0.00;

    k31 = 0.00;
    k32 = 0.00;
    k33 = 0.00;

    for(int i = 0;i < NProtein;i++) 
   {
    k11 += (double)kmat_11_host[i];
    k12 += (double)kmat_12_host[i];
    k13 += (double)kmat_13_host[i];

    k21 += (double)kmat_21_host[i];
    k22 += (double)kmat_22_host[i];
    k23 += (double)kmat_23_host[i];

    k31 += (double)kmat_31_host[i];
    k32 += (double)kmat_32_host[i];
    k33 += (double)kmat_33_host[i];

   }

   *one1 = k11;
   *one2 = k12;
   *one3 = k13;
   *two1 = k21;
   *two2 = k22;
   *two3 = k23;
   *three1 = k31;
   *three2 = k32;
   *three3 = k33;

}


/*

calc_kmat_sel_gpu_: (called from fortran as "calc_kmat_sel_gpu") is the fortran interface function for the 
                calculation of the 9 elements of the Kabsch matrix.  This routine assumes the selected crystal
                structure that we will test against has already been copied to the GPU.

*/
extern "C" void calc_kmat_sel_gpu_(float* xyz, double *one1, double *one2, double *one3, double *two1, double *two2, double *two3, double *three1, double *three2, double *three3)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock;
   if(NSelection >= 1000) {
     numThreadsPerBlock = 32;
   } else {
     numThreadsPerBlock = 256;
   }

   int numBlocks = NSelection / numThreadsPerBlock + (NSelection % numThreadsPerBlock == 0?0:1);

   float* xx;
   float* yy;
   float* zz;
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];

    int j = -3;
    for(int i = 0;i < NSelection;i++) 
   {
    j+=3;
    // fill up the host with the protein position at the current timestep
    prot_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);
   }

    // copy the current protein position to the device
   cudaMemcpy(prot_dev, prot_host, memSelectionF4, cudaMemcpyHostToDevice);


   // calculate the individual Kabsch elements on the GPU
   calckMat_GPU<<<numBlocks, numThreadsPerBlock>>>(NSelection,prot_dev,prot_crys_dev,kmat_11_dev, kmat_12_dev, kmat_13_dev, kmat_21_dev, kmat_22_dev, kmat_23_dev, kmat_31_dev, kmat_32_dev, kmat_33_dev);

   // copy the element arrays of the Kabsch matrix to the host
   cudaMemcpy(kmat_11_host, kmat_11_dev, memKabschSelF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_12_host, kmat_12_dev, memKabschSelF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_13_host, kmat_13_dev, memKabschSelF1, cudaMemcpyDeviceToHost);

   cudaMemcpy(kmat_21_host, kmat_21_dev, memKabschSelF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_22_host, kmat_22_dev, memKabschSelF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_23_host, kmat_23_dev, memKabschSelF1, cudaMemcpyDeviceToHost);

   cudaMemcpy(kmat_31_host, kmat_31_dev, memKabschSelF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_32_host, kmat_32_dev, memKabschSelF1, cudaMemcpyDeviceToHost);
   cudaMemcpy(kmat_33_host, kmat_33_dev, memKabschSelF1, cudaMemcpyDeviceToHost);

    double k11, k12, k13;
    double k21, k22, k23;
    double k31, k32, k33;

    k11 = 0.00;
    k12 = 0.00;
    k13 = 0.00;

    k21 = 0.00;
    k22 = 0.00;
    k23 = 0.00;

    k31 = 0.00;
    k32 = 0.00;
    k33 = 0.00;

    for(int i = 0;i < NSelection;i++) 
   {
    k11 += (double)kmat_11_host[i];
    k12 += (double)kmat_12_host[i];
    k13 += (double)kmat_13_host[i];

    k21 += (double)kmat_21_host[i];
    k22 += (double)kmat_22_host[i];
    k23 += (double)kmat_23_host[i];

    k31 += (double)kmat_31_host[i];
    k32 += (double)kmat_32_host[i];
    k33 += (double)kmat_33_host[i];

   }

   *one1 = k11;
   *one2 = k12;
   *one3 = k13;
   *two1 = k21;
   *two2 = k22;
   *two3 = k23;
   *three1 = k31;
   *three2 = k32;
   *three3 = k33;

}

/*

get_avg_sel_gpu_: (called from fortran as "accum_avg_sel_gpu") is the fortran interface function for the 
               routine to retrieve the average structure of the selection from the GPU

*/
extern "C" void get_avg_sel_gpu_(float *xavg,float *yavg,float *zavg)
{
   

   float* xx;
   float* yy;
   float* zz;
   xx = &xavg[0];
   yy = &yavg[0];
   zz = &zavg[0];

    // copy the current protein position to the host
   cudaMemcpy(prot_avg_host, prot_avg_dev, memSelectionF4, cudaMemcpyDeviceToHost);

    for(int i = 0;i < NSelection;i++) 
   {
    xx[i] = prot_avg_host[i].x;
    yy[i] = prot_avg_host[i].y;
    zz[i] = prot_avg_host[i].z;
   }


}


/*

get_avg_gpu_: (called from fortran as "accum_avg_gpu") is the fortran interface function for the 
               routine to retrieve the average structure from the GPU

*/
extern "C" void get_avg_gpu_(float *xavg,float *yavg,float *zavg)
{
   

   float* xx;
   float* yy;
   float* zz;
   xx = &xavg[0];
   yy = &yavg[0];
   zz = &zavg[0];

    // copy the current protein position to the host
   cudaMemcpy(prot_avg_host, prot_avg_dev, memProteinF4, cudaMemcpyDeviceToHost);

    for(int i = 0;i < NProtein;i++) 
   {
    xx[i] = prot_avg_host[i].x;
    yy[i] = prot_avg_host[i].y;
    zz[i] = prot_avg_host[i].z;
   }


}

/*

get_covar_gpu_: (called from fortran as "accum_avg_gpu") is the fortran interface function for the 
               routine to retrieve the covariance matrix from the GPU

*/
extern "C" void get_covar_gpu_(float *covar_cpu)
{
   
    // copy the covariance matrix directory to fortrans memory space..... 
    // ... 
    // ...
    // this might be a mistake
   cudaMemcpy(covar_cpu, covar_dev, memProteinDOF * NProteinDOF, cudaMemcpyDeviceToHost);

   /// QUESTION: DO WE NEED TO TRANSPOSE THIS PUPPY???
   /// ANSWER: NO, ITS SYMMETRIC ABOUT THE DIAGONAL

}

/*

get_covar_sel_gpu_: (called from fortran as "get_covar_sel_gpu") is the fortran interface function for the 
               routine to retrieve the covariance matrix for selected atoms from the GPU

*/
extern "C" void get_covar_sel_gpu_(float *covar_cpu)
{
   
    // copy the covariance matrix directory to fortrans memory space..... 
    // ... 
    // ...
    // this might be a mistake
   cudaMemcpy(covar_cpu, covar_dev, memSelectionDOF * NSelectionDOF, cudaMemcpyDeviceToHost);

   /// QUESTION: DO WE NEED TO TRANSPOSE THIS PUPPY???
   /// ANSWER: NO, ITS SYMMETRIC ABOUT THE DIAGONAL

}

/*

store_crys_gpu_: (called from fortran as "store_crys_gpu") is the fortran interface function for the 
               routine to store the crystal structure along with the masses on the GPU

*/
extern "C" void store_crys_gpu_(float *xyz, float* qrms,bool* setVectors)
{

   float* xx;
   float* yy;
   float* zz;
   float* protMass;

   protMass = &qrms[2+iPStart]; /// points to 2 for the mass position
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];
    
    int j = -3;
    int k = -1;
    float sqrtM;
    for(int i = 0;i < NProtein;i++) 
   {
    j+=3;
    prot_crys_host[i].x = xx[j];
    prot_crys_host[i].y = yy[j];
    prot_crys_host[i].z = zz[j];
    prot_crys_host[i].w = protMass[j];

    if(*setVectors)
    {
     sqrtM = sqrtf(protMass[j]);
     k++;
     prot_avg_vec_host[k] = xx[j]*sqrtM;     
     k++;
     prot_avg_vec_host[k] = yy[j]*sqrtM;
     k++;
     prot_avg_vec_host[k] = zz[j]*sqrtM;
    }

//     if(i == (NProtein - 1))
//     {
//       printf("HOST(CRYS INFO)  \n");
//       printf("sqrtM = %f,  \n",sqrtM);
//       printf("prot_avg_vec_host[k,k-1,k-2] = %f, %f, %f \n",prot_avg_vec_host[k],prot_avg_vec_host[k-1],prot_avg_vec_host[k-2]); 
//       printf("prot_crys_host[i].x,y,z = %f, %f, %f \n",prot_crys_host[i].x,prot_crys_host[i].y,prot_crys_host[i].z); 
//     } 

   }

   // copy the crystal structure position is copied to the device
   cudaMemcpy(prot_crys_dev, prot_crys_host, memProteinF4, cudaMemcpyHostToDevice);

    if(*setVectors)
    {
     // copy the mass weighted average position is copied to the device in vector format 
     cudaMemcpy(prot_avg_vec_dev, prot_avg_vec_host, memProteinDOF, cudaMemcpyHostToDevice);

    }
}

/*

store_crys_sel_gpu_: (called from fortran as "store_crys_sel_gpu") is the fortran interface function for the 
               routine to store the selected crystal structure along with the masses on the GPU

*/
extern "C" void store_crys_sel_gpu_(float *xyz, float* qrms,bool* setVectors,int* selIndexes)
{

   float* xx;
   float* yy;
   float* zz;
   float* protMass;

   protMass = &qrms[2]; /// points to 2 for the mass position
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];
    
    int j,k;
    int l = -1;
    int m = -3;
    float sqrtM;
    for(int i = 0;i < NSelection;i++) 
   {
    k = selIndexes[i];
    j = 3*k - 3;

    m+=3;
    prot_crys_host[i].x = xx[m];
    prot_crys_host[i].y = yy[m];
    prot_crys_host[i].z = zz[m];
    prot_crys_host[i].w = protMass[j];

    if(*setVectors)
    {
     sqrtM = sqrtf(protMass[j]);
     l++;
     prot_avg_vec_host[l] = xx[m]*sqrtM;     
     l++;
     prot_avg_vec_host[l] = yy[m]*sqrtM;
     l++;
     prot_avg_vec_host[l] = zz[m]*sqrtM;
    }

//     if(i == (NProtein - 1))
//     {
//       printf("HOST(CRYS INFO)  \n");
//       printf("sqrtM = %f,  \n",sqrtM);
//       printf("prot_avg_vec_host[k,k-1,k-2] = %f, %f, %f \n",prot_avg_vec_host[k],prot_avg_vec_host[k-1],prot_avg_vec_host[k-2]); 
//       printf("prot_crys_host[i].x,y,z = %f, %f, %f \n",prot_crys_host[i].x,prot_crys_host[i].y,prot_crys_host[i].z); 
//     } 

   }

   // copy the crystal structure position is copied to the device
   cudaMemcpy(prot_crys_dev, prot_crys_host, memSelectionF4, cudaMemcpyHostToDevice);

    if(*setVectors)
    {
     // copy the mass weighted average position is copied to the device in vector format 
     cudaMemcpy(prot_avg_vec_dev, prot_avg_vec_host, memSelectionDOF, cudaMemcpyHostToDevice);

    }
}


/*

get_crys_gpu_: (called from fortran as "get_crys_gpu") is the fortran interface function for the 
               routine to pull the crystal structure from the GPU to the CPU

*/
extern "C" void get_crys_gpu_(float *xyz)
{

   float* xx;
   float* yy;
   float* zz;
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];

   // copy the current protein position to the host
   cudaMemcpy(prot_crys_host, prot_crys_dev, memProteinF4, cudaMemcpyDeviceToHost);

    int j = -3;
    for(int i = 0;i < NProtein;i++) 
   {
    j+=3;
    xx[j] = prot_crys_host[i].x;
    yy[j] = prot_crys_host[i].y;
    zz[j] = prot_crys_host[i].z;
   }

}

void checkCUDAError22(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}

