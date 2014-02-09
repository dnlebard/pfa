#include <stdio.h>
#include <stdbool.h>
#include <cuda.h>

/*********************************
HEADER SECTION
//////////////////////////////////
*********************************/
//! Texture for reading water parameters
//texture<float4, 1, cudaReadModeElementType> water_tex;

//! Pointers for storing the water and energy arrays
float4* water_dev; // water on the device
float4* water_host;   // float4 version of the water on the host

// System-wide number of waters
int NWaterAtoms;
int iWatNum; // the starting atomic index
int iWatEnd; // the ending atomic index

// Size infos
//size_t memFloat4;
//size_t memFloat;
unsigned int memFloat4;
unsigned int memFloat;

float* water_en_dev; // water energy on the device
float* water_en_host; // water energy on the device

//////////////////////////////////
// Simple utility function to check for CUDA runtime errors
void checkCUDAError(const char* msg);

extern "C" {extern "C" void alloc_gpu_dvos_(int *rank, int* iWatStart, int *iAtoms,int* iGPUs,bool *doGPUInit);}

extern "C" { extern "C" void dealloc_gpu_dvos_();}

extern "C" { extern "C" void dvos_gpu_(float* pX, float* pY, float* pZ, float *dQ, float *boxX, float *boxY, float *boxZ, float *xyz, float* qrms, float* dvos_i); 
}
/*********************************
END HEADER SECTION
//////////////////////////////////
*********************************/
__global__ void dVos_i_GPU(float pX, float pY, float pZ, float dQ, float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int N, bool isTruncated, float4* water_dev, float* water_en_dev)
{

     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < N)
     {// then we actually do something...


   // get our water's atom's position
    float4 water_pos = water_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   

   float wX = water_pos.x;
   float wY = water_pos.y;
   float wZ = water_pos.z;
   float wQ = water_pos.w; // charge

   float dX =  pX - wX;
   dX -= boxX * rintf((dX * boxXinv));

   float dY =  pY - wY;
   dY -= boxY * rintf((dY * boxYinv));

   float dZ =  pZ - wZ;
   dZ -= boxZ * rintf((dZ * boxZinv));

   float rWP2 = dX*dX + dY*dY + dZ*dZ;

   float inv_rWP = rsqrtf(rWP2);

   float dvos_i = inv_rWP * wQ * dQ;

   water_en_dev[idx_water] = dvos_i; // update the energy array

   }// else we just sync our threads and go...

   __syncthreads();

     return; 

}

extern "C" void alloc_gpu_dvos_(int *rank, int* iWatStart, int *iAtoms,int* iGPUs,bool *doGPUInit)
{

   iWatNum = (*iWatStart * 3) - 3; // the starting atomic index
   iWatEnd = (*iAtoms - 1) * 3; // the ending atomic index
   NWaterAtoms = *iAtoms - *iWatStart + 1;
   // Find our local GPU rank on the node
    char cudaOut[100];
    int iChr;
    int rank_gpu = *rank % *iGPUs;

   if(*doGPUInit)
   {


    //printf("===== SETTING THE CUDA DEVICE ON RANK= %i to CUDA_RANK= %i \n",*rank,rank_gpu);

    cudaSetDevice(rank_gpu); 
    //printf("====== DONE === SETTING THE CUDA DEVICE ON RANK= %i to CUDA_RANK= %i \n",*rank,rank_gpu);

    iChr = sprintf(cudaOut,"cudaSetDevice -- setting device to %i on rank %i \n",rank_gpu,*rank);
    if(iChr < 0 ){printf("SHOULD I CARE?\n");}
    checkCUDAError(cudaOut);
    *doGPUInit = false;
   }

   //char *methodName;
   //methodName = "alloc_gpu_dvos";
   //alloc_gpu_rank_(rank, iGPUs, methodName);

   // Find our local GPU rank on the node
   //int rank_gpu = *rank % *iGPUs;
   //cudaSetDevice(rank_gpu); // THIS COULD BE A DANGEROUS PLACE TO MAKE THIS CALL

  //cudaSetDevice(*rank);
  // checkCUDAError("cudaSetDevice -- setting device to 0");

    //printf("===== SETTING THE CUDA SIZES ON RANK= %i to CUDA_RANK= %i \n",*rank,rank_gpu);
   //------------------------------------------
   // PUT THIS IN AN ALLOCATION ROUTINE
   //------------------------------------------
   memFloat4 = sizeof(float4) *  NWaterAtoms ; 
   memFloat = sizeof(float) *  NWaterAtoms ;
   // printf("=========>. FOUND SIZES RANK= %i to CUDA_RANK= %i \n",*rank,rank_gpu);
   water_host = (float4*)malloc(memFloat4);
   // printf("=========>. MALLOCED WATER HOST = %i to CUDA_RANK= %i \n",*rank,rank_gpu);
   water_en_host = (float*)malloc(memFloat);
  //  printf("=========>. MALLOCED WATER EN HOST = %i to CUDA_RANK= %i \n",*rank,rank_gpu);

   cudaMalloc( (void **) &water_dev, memFloat4);
   checkCUDAError("cudaMalloc -- water_dev");
  //  printf("=========>. CUDA MALLOCED WATER DEV = %i to CUDA_RANK= %i \n",*rank,rank_gpu);

   cudaMalloc( (void **) &water_en_dev,memFloat);
   checkCUDAError("cudaMalloc -- water_en_dev");
  //  printf("=========>. CUDA MALLOCED WATER EN DEV = %i to CUDA_RANK= %i \n",*rank,rank_gpu);

  //  printf("===== DONE === SETTING THE CUDA SIZES ON RANK= %i to CUDA_RANK= %i \n",*rank,rank_gpu);

   //------------------------------------------
   // END ALLOCATION ROUTINE
   //------------------------------------------

}

extern "C" void dealloc_gpu_dvos_()
{
   //------------------------------------------
   // PUT THIS IN A DEALLOCATION ROUTINE
   //------------------------------------------
   // free device memory
    cudaFree(water_dev);
    cudaFree(water_en_dev);

   // free up host memory
      delete[] water_host;
      water_host = NULL;

      delete[] water_en_host;
      water_en_host = NULL;

   //------------------------------------------
   // END DEALLOCATION ROUTINE
   //------------------------------------------

}

/*
myc_ : a code to call the cuda kernel from fortran
*/
extern "C" void dvos_gpu_(float* pX, float* pY, float* pZ, float *dQ, float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvos_i)
{


   bool isTruncated = false;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256/4;
   int numBlocks = NWaterAtoms / numThreadsPerBlock + (NWaterAtoms % numThreadsPerBlock == 0?0:1);


   float* xx;
   float* yy;
   float* zz;
   float* qq;
   float boxXinv, boxYinv, boxZinv;

 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);

/*
   xx = &xyz[(*iWatStart * 3) -3];
   yy = &xyz[(*iWatStart * 3) -2];
   zz = &xyz[(*iWatStart * 3) -1];
   qq = &qrms[(*iWatStart * 3) -3];
*/
   xx = &xyz[iWatNum];
   yy = &xyz[iWatNum+1];
   zz = &xyz[iWatNum+2];
   qq = &qrms[iWatNum];

   int j = -3;
   // loop over all waters and fill up the water array
   for(int i = 0;i < NWaterAtoms;i++) // should this just be NWaterAtoms????
   {
    j+=3;
    // fill up the water on the host in float4s
    water_host[i] = make_float4(xx[j],yy[j],zz[j],qq[j]);

   }

    // copy the water from the host to the device
   cudaMemcpy(water_dev, water_host, memFloat4, cudaMemcpyHostToDevice);
   //checkCUDAError("cudaMemcpy -- host to device");

   dVos_i_GPU<<<numBlocks, numThreadsPerBlock>>>(*pX,*pY,*pZ,*dQ,
                                   *boxX,*boxY,*boxZ,
                                   boxXinv,boxYinv,boxZinv,
                                   NWaterAtoms,isTruncated,
                                   water_dev,water_en_dev);

  // copy the water from the host to the device
   cudaMemcpy(water_en_host,water_en_dev,memFloat,cudaMemcpyDeviceToHost);
   //checkCUDAError("cudaMemcpy -- results from device to host");


  // Sum the whole array
    float dVos_i_sum = 0.0f;
    for(int i=0;i<NWaterAtoms;i++)
     dVos_i_sum += water_en_host[i];

    *dvos_i = dVos_i_sum;

}


void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}

