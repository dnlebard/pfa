#include <stdio.h>
#include <cuda.h>

/*********************************
HEADER SECTION
//////////////////////////////////
*********************************/
//! Texture for reading water parameters
//texture<float4, 1, cudaReadModeElementType> water_tex;

//! Pointers for storing the water, protein, and shell
float4* water_oxy_dev;    // water's oxygen on the device
float4* water_oxy_host;   // water's oxygen on the host

float4* water_h1_dev;    // water's hydrogen #1 on the device
float4* water_h1_host;   // water's hydrogen #1 on the host

float4* water_h2_dev;    // water's hydrogen #2 on the device
float4* water_h2_host;   // water's hydrogen #2 on the host

float4* water_binfo_host;    // the bin-cos float4 storage (note: z,w are empty, useable) 
float4* water_binfo_dev;    // the bin-cos float4 storage (note: z,w are empty, useable) 

float4* protein_dev;  // protein atoms on the device
float4* protein_host; // protein atoms on the host

float4* protein_sel_dev;  // selected protein atoms on the device
float4* protein_sel_host; // selected protein atoms on the host

float4* dvos_dev;  // dV0s atoms on the device
float4* dvos_host; // dV0s atoms on the host

float4* protMu_dev;  // protein dipole on the device
float4* protMu_host; // protein dipole on the host
                      //
float4* protMu_sel_dev;  // selected protein dipole on the device
float4* protMu_sel_host; // selected protein dipole on the host
                     //
float4* mc_ran_host;  // the random numbers for r(x,y,z) generation
float4* mc_ran_dev;   // the random numbers for r(x,y,z) generation
                      //
                      //
float* mc_prot_host;   // (output) counts of test positions that fall in the protein volume
float* mc_prot_dev;   // (output) counts of test positions that fall in the protein volume
float* mc_wat_host;    // (output) counts of test positions that fall in the water volume
float* mc_wat_dev;    // (output) counts of test positions that fall in the water volume
float* mc_count_host;  // (output) counts of test positions that fall in the sphere
float* mc_count_dev;  // (output) counts of test positions that fall in the sphere
                       //
float* p2arg_host;    // the hosts array of p2arg's for all 1st shell waters
float* p2arg_dev;     // the device array of p2arg's for all 1st shell waters
                      //
float* dVFull_host;    // the hosts array of dV for all 1st shell waters/all proteins
float* dVFull_dev;     // the device array of dV for all 1st shell waters/all proteins
                       //
float* dVDelta_host;    // the hosts array of dV0s for all 1st shell waters/all proteins
float* dVDelta_dev;     // the device array of dV0s for all 1st shell waters/all proteins

float* dVpol_host;      // the hosts array of the Vpol correction
float* dVpol_dev;      // the device array of the Vpol correction

float* dQpol_host;      // the hosts array of the Qpol correction
float* dQpol_dev;      // the device array of the Qpol correction
                        //
                       //
int*  lshell_dev;     // The lshell array on the device
int* lshell_host;     // the lshell array on the host (long and unsorted)
                      //
// NOTE: E-field is calculated for at water's oxygen only
float* water_ef_dev;    // water's electric field on the device
float* water_ef_host;   // water's electric field on the host
////////////////////////

// System-wide numbers
int NwAtoms; // the number of water atoms
int NDVAtoms; // the number of DV0s atoms
int NMCAttempts; // the number of attempts on the mc volume calc
int NProtAtoms; // the number of protein atoms
int NSelAtoms; // the number of selected atoms
int NShellMax; // the maximum number of shell atoms
int iWaterIdx; // the starting water atomic index
int iLastWat; // the ending atomic index
int iFortWatIdx; // the fortran water starting index
int iProtEnd; // the ending protein index
int iDVEnd; // the ending DVAtom index

// Size infos
unsigned int memWaterF4;
unsigned int memProtF4;
unsigned int memSelF4;
unsigned int memFloatWat;
unsigned int memShellInt;
unsigned int memDVAtomsF4;
unsigned int memMCVolF1;
unsigned int memMCVolF4;

//////////////////////////////////
// Simple utility function to check for CUDA runtime errors
void checkCUDAError2(const char* msg);

extern "C" { extern "C" void alloc_gpu_shell_(int* rank, int* iWatStart, int* iPrEnd,int* iAtoms,int* iDVAtoms, int* iAttempts, int* iSelect, int* iGPUs, bool* doGPUInit);}

extern "C" { extern "C" void dealloc_gpu_shell_();}

extern "C" { extern "C" void get_shell_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* rCut,int* lShell,int* nShell, bool* checkForAll);}

extern "C" { extern "C" void get_p2_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms,int* lShell,int* nShell,float* p2arg, float* p2Shell);}

extern "C" {extern "C" extern "C" void get_prmu_gpu_(float* qrms, float* ddX,float* ddY, float* ddZ, 
                         float* prX,float* prY, float* prZ, float* totMorQ, int* massOrQ, bool* doOutPos, float* outputPos);}

extern "C" {extern "C" void get_prmu_sel_gpu_(float* qrms, float* ddX,float* ddY, float* ddZ, 
                              float* prX,float* prY, float* prZ, float* totMorQ, int* massOrQ, int* selIndices, bool* doOutPos, float* outputPos);}

extern "C" {extern "C" void get_shl_cos_gpu_(float *xyz, float* qrms, float* prMx, float* prMy, float* prMz, float* dCos, int* lShell,int* nShell,float* cosShell,int* binCos);}

extern "C" {extern "C" void get_dv0s_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvAtoms, float* dvQs, int* lShell,int* nShell,float* dVFull, float* dVDelta);}

extern "C" {extern "C" void get_efld_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvAtoms, float* dvQs, int* lShell,int* nShell,float* efield);}

extern "C" {extern "C" void get_vpol_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvAtoms, float* dvQs, int* lShell,int* nShell,float* dVpol);}

extern "C" {extern "C" void get_qpol_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, int* lShell,int* nShell,float* dQpol);}

extern "C" {extern "C" void get_vol_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, int* lShell,int* nShell, float* mcRadius, float* randx, float* randy, float* randz, float* centerX, float* centerY, float* centerZ, float* volProt, float* volWat, float* epsProt, float* epsWat);}

extern "C" {extern "C" void set_prot_sel_gpu_(float* xyz);}

extern "C" {extern "C" void set_prot_gpu_(float *xyz);}
/*********************************
END HEADER SECTION
//////////////////////////////////
*********************************/
/*
   qpol_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  dQpol_dev);
*/
__global__ void qpol_shl_GPU(float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int nSh, bool isTruncated, float4* water_oxy_dev, float4* water_h1_dev, float4* water_h2_dev, float* dQpol_dev)
{
     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < nSh)
     {// then we actually do something...
  
      // get our oxygen's atom's position
      float4 water_o_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float oX = water_o_pos.x;
      float oY = water_o_pos.y;
      float oZ = water_o_pos.z;
      float oQ = water_o_pos.w; 

   // get our hydrogen's #1 atom's position
      float4 water_h1_pos = water_h1_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h1X = water_h1_pos.x;
      float h1Y = water_h1_pos.y;
      float h1Z = water_h1_pos.z;
      float h1Q = water_h1_pos.w; 

   // get our hydrogen's #1 atom's position
      float4 water_h2_pos = water_h2_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h2X = water_h2_pos.x;
      float h2Y = water_h2_pos.y;
      float h2Z = water_h2_pos.z;
      float h2Q = water_h2_pos.w; 

       // oxygen
       float dXo =  oX;
       dXo -= boxX * rintf((dXo * boxXinv));
       float dYo =  oY;
       dYo -= boxY * rintf((dYo * boxYinv));
       float dZo =  oZ;
       dZo -= boxZ * rintf((dZo * boxZinv));
       //float rWoP2 = dXo*dXo + dYo*dYo + dZo*dZo;
       //float rWoP_inv = rsqrtf(rWoP2);
       //float enWo = oQ*pQ*rWoP2;

       // h1
       float dXh1 =  h1X;
       dXh1 -= boxX * rintf((dXh1 * boxXinv));
       float dYh1 =  h1Y;
       dYh1 -= boxY * rintf((dYh1 * boxYinv));
       float dZh1 =  h1Z;
       dZh1 -= boxZ * rintf((dZh1 * boxZinv));

       float dXoh1 = dXo - dXh1;
       float dYoh1 = dYo - dYh1;
       float dZoh1 = dZo - dZh1;
       float roh12 = dXoh1*dXoh1 + dYoh1*dYoh1 + dZoh1*dZoh1;
       float roh12_inv = rsqrtf(roh12);
       //float rWh1P2 = dXh1*dXh1 + dYh1*dYh1 + dZh1*dZh1;
       //float rWh1P_inv = rsqrtf(rWh1P2);
       //float enWh1 = h1Q*pQ*rWh1P2;

       // h2
       float dXh2 =  h2X;
       dXh2 -= boxX * rintf((dXh2 * boxXinv));
       float dYh2 =  h2Y;
       dYh2 -= boxY * rintf((dYh2 * boxYinv));
       float dZh2 =  h2Z;
       dZh2 -= boxZ * rintf((dZh2 * boxZinv));

       float dXoh2 = dXo - dXh2;
       float dYoh2 = dYo - dYh2;
       float dZoh2 = dZo - dZh2;
       float roh22 = dXoh2*dXoh2 + dYoh2*dYoh2 + dZoh2*dZoh2;
       float roh22_inv = rsqrtf(roh22);
       //float rWh2P2 = dXh2*dXh2 + dYh2*dYh2 + dZh2*dZh2;
       //float rWh2P_inv = rsqrtf(rWh2P2);
       //float enWh2 = h2Q*pQ*rWh2P2;

      unsigned int keepNew;  
      keepNew = 1;
      if( roh12_inv > 0.40000f || roh22_inv > 0.40000f) keepNew = 0;

      float woX;
      float woY;
      float woZ;

      float wmX;
      float wmY;
      float wmZ;
      wmX = 0.0f;
      wmY = 0.0f;
      wmZ = 0.0f;

      
      if(keepNew == 1) 
      {
       woX = dXo;
       woY = dYo;
       woZ = dZo;

       wmX = oQ*dXo;
       wmY = oQ*dYo;
       wmZ = oQ*dZo;

       wmX += h1Q*dXh1;
       wmY += h1Q*dYh1;
       wmZ += h1Q*dZh1;

       wmX += h2Q*dXh2;
       wmY += h2Q*dYh2;
       wmZ += h2Q*dZh2;

      }else {
       woX = oX;
       woY = oY;
       woZ = oZ;

       wmX = oQ*oX;
       wmY = oQ*oY;
       wmZ = oQ*oZ;

       wmX += h1Q*h1X;
       wmY += h1Q*h1Y;
       wmZ += h1Q*h1Z;

       wmX += h2Q*h2X;
       wmY += h2Q*h2Y;
       wmZ += h2Q*h2Z;
      }


       float RoMw = woX*wmX + woY*wmY + woZ*wmZ;
       dQpol_dev[idx_water] = RoMw;

     }
   __syncthreads();
     return; 



}

/*
   vpol_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NDVAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  dvos_dev,
                                                  dVpol_dev);
*/
__global__ void vpol_shl_GPU(float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int nSh, int NdV, bool isTruncated, float4* water_oxy_dev, float4* water_h1_dev, float4* water_h2_dev, float4* dvos_dev, float* dVpol_dev)
{
     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < nSh)
     {// then we actually do something...
  
      // get our oxygen's atom's position
      float4 water_o_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float oX = water_o_pos.x;
      float oY = water_o_pos.y;
      float oZ = water_o_pos.z;
      float oQ = water_o_pos.w; 

   // get our hydrogen's #1 atom's position
      float4 water_h1_pos = water_h1_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h1X = water_h1_pos.x;
      float h1Y = water_h1_pos.y;
      float h1Z = water_h1_pos.z;
      float h1Q = water_h1_pos.w; 

   // get our hydrogen's #1 atom's position
      float4 water_h2_pos = water_h2_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h2X = water_h2_pos.x;
      float h2Y = water_h2_pos.y;
      float h2Z = water_h2_pos.z;
      float h2Q = water_h2_pos.w; 

      float en_dV = 0.0f;
      for(int i=0;i<NdV;i++)
      {
       float4 dvos_pos = dvos_dev[i];
       float pX = dvos_pos.x;
       float pY = dvos_pos.y;
       float pZ = dvos_pos.z; 
       float pQ = dvos_pos.w;

       // oxygen
       float dXo =  pX - oX;
       dXo -= boxX * rintf((dXo * boxXinv));
       float dYo =  pY - oY;
       dYo -= boxY * rintf((dYo * boxYinv));
       float dZo =  pZ - oZ;
       dZo -= boxZ * rintf((dZo * boxZinv));
       float rWoP2 = dXo*dXo + dYo*dYo + dZo*dZo;
       //float rWoP_inv = rsqrtf(rWoP2);
       float enWo = oQ*pQ*rWoP2;

       // h1
       float dXh1 =  pX - h1X;
       dXh1 -= boxX * rintf((dXh1 * boxXinv));
       float dYh1 =  pY - h1Y;
       dYh1 -= boxY * rintf((dYh1 * boxYinv));
       float dZh1 =  pZ - h1Z;
       dZh1 -= boxZ * rintf((dZh1 * boxZinv));
       float rWh1P2 = dXh1*dXh1 + dYh1*dYh1 + dZh1*dZh1;
       //float rWh1P_inv = rsqrtf(rWh1P2);
       float enWh1 = h1Q*pQ*rWh1P2;

       // h2
       float dXh2 =  pX - h2X;
       dXh2 -= boxX * rintf((dXh2 * boxXinv));
       float dYh2 =  pY - h2Y;
       dYh2 -= boxY * rintf((dYh2 * boxYinv));
       float dZh2 =  pZ - h2Z;
       dZh2 -= boxZ * rintf((dZh2 * boxZinv));
       float rWh2P2 = dXh2*dXh2 + dYh2*dYh2 + dZh2*dZh2;
       //float rWh2P_inv = rsqrtf(rWh2P2);
       float enWh2 = h2Q*pQ*rWh2P2;

       en_dV+=(enWo + enWh1 + enWh2);

     } // end for

       dVpol_dev[idx_water] = en_dV;

     }
   __syncthreads();
     return; 



}

/*

      efield_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NDVAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,dvos_dev,
                                                  water_ef_dev);

*/
__global__ void efield_shl_GPU(float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int nSh, int NdV, bool isTruncated, float4* water_oxy_dev, float4* dvos_dev, float* water_ef_dev)
{
     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < nSh)
     {// then we actually do something...
  
      // get our oxygen's atom's position
      float4 water_o_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float oX = water_o_pos.x;
      float oY = water_o_pos.y;
      float oZ = water_o_pos.z;


      float efx_dV = 0.0f;
      float efy_dV = 0.0f;
      float efz_dV = 0.0f;
      for(int i=0;i<NdV;i++)
      {
       float4 dvos_pos = dvos_dev[i];
       float pX = dvos_pos.x;
       float pY = dvos_pos.y;
       float pZ = dvos_pos.z; 
       float pQ = dvos_pos.w;

       // oxygen
       float dXo =  oX - pX;
       dXo -= boxX * rintf((dXo * boxXinv));
       float dYo =  oY - pY;
       dYo -= boxY * rintf((dYo * boxYinv));
       float dZo =  oZ - pZ;
       dZo -= boxZ * rintf((dZo * boxZinv));
       float rWoP2 = dXo*dXo + dYo*dYo + dZo*dZo;
       float rWoP_inv = rsqrtf(rWoP2);
       float rWoP_inv3 = rWoP_inv*rWoP_inv*rWoP_inv;

       float r3dq = pQ*rWoP_inv3;

       float efx = dXo*r3dq;
       float efy = dYo*r3dq;
       float efz = dZo*r3dq;

       efx_dV+=efx;
       efy_dV+=efy;
       efz_dV+=efz;

     } // end for

       float efx2 = efx_dV*efx_dV;
       float efy2 = efy_dV*efy_dV;
       float efz2 = efz_dV*efz_dV;

       float efMag2 = efx2+efy2+efz2;
       float efMag = sqrtf(efMag2);

       water_ef_dev[idx_water] = efMag;

     }
     __syncthreads();
     return; 



}

/*
   
   vol_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  *mcRadius,
                                                  *centerX, *centerY, *centerZ,
                                                  nSh,NProtAtoms,NMCAttempts,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,water_h2_dev,
                                                  protein_dev,
                                                  mc_ran_dev,
                                                  mc_prot_dev,mc_wat_dev,mc_count_dev);


*/
__global__ void vol_shl_GPU(float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, float totRadius, float ctrX, float ctrY, float ctrZ, int nSh, int Np, int NMoves, bool isTruncated, float4* water_oxy_dev, float4* water_h1_dev, float4* water_h2_dev, float4* protein_dev, float4* mc_ran_dev, float* mc_prot_dev, float* mc_wat_dev, float* mc_count_dev)
{
     unsigned int idx_mc = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_mc < NMoves)
     {// then we actually do something...

      // get our test position
      float4 ran_num = mc_ran_dev[idx_mc];
      float rX = ran_num.x;
      float rY = ran_num.y;
      float rZ = ran_num.z;

      float x = (2.0f*rX - 1.0f)*totRadius;
      float y = (2.0f*rY - 1.0f)*totRadius;
      float z = (2.0f*rZ - 1.0f)*totRadius;

      float xyz2 = x*x + y*y + z*z;

      float totR2 = totRadius*totRadius;


      float countVal;
      float protVal= 0.0f;


      if (xyz2 < totR2)
      { 
       countVal = 1.0f;


/*
        printf("x,y,z = %f,%f,%f \n",x,y,z);
        printf("cx,cy,cz = %f,%f,%f \n",ctrX,ctrY,ctrZ);
        printf("xyz2 = %f \n",xyz2);
        printf("totR2 = %f \n",totR2);
        printf("countVal = %f \n",countVal);
*/
      } else {
       countVal = 0.0f;
      }


      for(int i=0;i<Np;i++)
      {
       float4 prot_pos = protein_dev[i];
       float pX = prot_pos.x;
       float pY = prot_pos.y;
       float pZ = prot_pos.z; 
       float pR = prot_pos.w;
       float r2 = pR * pR;

       // check if we fall inside our protein
       float dpX =  pX - ctrX;
       dpX -= x;
       float dpY =  pY - ctrY;
       dpY -= y;
       float dpZ =  pZ - ctrZ;
       dpZ -= z;
       float dp2 = dpX*dpX + dpY*dpY + dpZ*dpZ;

       if (dp2 < r2)
        protVal += 1.0f;

      } // end for


      float oVal = 0.0f;
      float h1Val = 0.0f;
      float h2Val = 0.0f;
      for(int i=0;i<nSh;i++)
      {
       float4 oxy_pos = water_oxy_dev[i];
       float oX = oxy_pos.x;
       float oY = oxy_pos.y;
       float oZ = oxy_pos.z; 
       float oR = oxy_pos.w;
       float r2 = oR * oR;
 
       // check if we fall inside our oxygen
       float doX =  oX - ctrX;
       doX -= boxX * rintf((doX * boxXinv));
       doX -= x;
       float doY =  oY - ctrY;
       doY -= boxY * rintf((doY * boxYinv));
       doY -= y;
       float doZ =  oZ - ctrZ;
       doZ -= boxZ * rintf((doZ * boxZinv));
       doZ -= z;
       float do2 = doX*doX + doY*doY + doZ*doZ;

       if(do2 < r2) 
       { oVal += 1.0f;}

       float4 h1_pos = water_h1_dev[i];
       float h1X = h1_pos.x;
       float h1Y = h1_pos.y;
       float h1Z = h1_pos.z; 
       float h1R = h1_pos.w;
       r2 = h1R * h1R;
 
       // check if we fall inside our first hydrogen
       float dh1X =  h1X - ctrX;
       dh1X -= boxX * rintf((dh1X * boxXinv));
       dh1X -= x;
       float dh1Y =  h1Y - ctrY;
       dh1Y -= boxY * rintf((dh1Y * boxYinv));
       dh1Y -= y;
       float dh1Z =  h1Z - ctrZ;
       dh1Z -= boxZ * rintf((dh1Z * boxZinv));
       dh1Z -= z;
       float dh12 = dh1X*dh1X + dh1Y*dh1Y + dh1Z*dh1Z;

       if(dh12 < r2) 
        { h1Val += 1.0f;}
  
       float4 h2_pos = water_h2_dev[i];
       float h2X = h2_pos.x;
       float h2Y = h2_pos.y;
       float h2Z = h2_pos.z; 
       float h2R = h2_pos.w;
       r2 = h2R * h2R;

       // check if we fall inside our second hydrogen
       float dh2X =  h2X - ctrX;
       dh2X -= boxX * rintf((dh2X * boxXinv));
       dh2X -= x;
       float dh2Y =  h2Y - ctrY;
       dh2Y -= boxY * rintf((dh2Y * boxYinv));
       dh2Y -= y;
       float dh2Z =  h2Z - ctrZ;
       dh2Z -= boxZ * rintf((dh2Z * boxZinv));
       dh2Z -= z;
       float dh22 = dh2X*dh2X + dh2Y*dh2Y + dh2Z*dh2Z;

       if(dh22 < r2) 
        { h2Val += 1.0f;}

      } // end for


      float watVal = oVal + h1Val + h2Val;        

      //if(protVal > 1.0) protVal = 1.0;
      if(watVal > 1.0) watVal = 1.0;


      mc_prot_dev[idx_mc] = protVal;
      mc_wat_dev[idx_mc] = watVal;
      mc_count_dev[idx_mc] = countVal;

     }
   __syncthreads();
     return; 


}

/*
   dvos_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NProtAtoms,NDVAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  protein_dev,dvos_dev,
                                                  dVFull_dev,dVDelta_dev);

*/
__global__ void dvos_shl_GPU(float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int nSh, int Np, int NdV, bool isTruncated, float4* water_oxy_dev, float4* water_h1_dev, float4* water_h2_dev, float4* protein_dev, float4* dvos_dev, float* dVFull_dev,float* dVDelta_dev)
{
     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < nSh)
     {// then we actually do something...
  
      // get our oxygen's atom's position
      float4 water_o_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float oX = water_o_pos.x;
      float oY = water_o_pos.y;
      float oZ = water_o_pos.z;
      float oQ = water_o_pos.w; 

   // get our hydrogen's #1 atom's position
      float4 water_h1_pos = water_h1_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h1X = water_h1_pos.x;
      float h1Y = water_h1_pos.y;
      float h1Z = water_h1_pos.z;
      float h1Q = water_h1_pos.w; 

   // get our hydrogen's #1 atom's position
      float4 water_h2_pos = water_h2_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h2X = water_h2_pos.x;
      float h2Y = water_h2_pos.y;
      float h2Z = water_h2_pos.z;
      float h2Q = water_h2_pos.w; 


      float en_Prot = 0.0f;
      for(int i=0;i<Np;i++)
      {
       float4 prot_pos = protein_dev[i];
       float pX = prot_pos.x;
       float pY = prot_pos.y;
       float pZ = prot_pos.z; 
       float pQ = prot_pos.w;

       // oxygen
       float dXo =  pX - oX;
       dXo -= boxX * rintf((dXo * boxXinv));
       float dYo =  pY - oY;
       dYo -= boxY * rintf((dYo * boxYinv));
       float dZo =  pZ - oZ;
       dZo -= boxZ * rintf((dZo * boxZinv));
       float rWoP2 = dXo*dXo + dYo*dYo + dZo*dZo;
       float rWoP_inv = rsqrtf(rWoP2);
       float enWo = oQ*pQ*rWoP_inv;

       // h1
       float dXh1 =  pX - h1X;
       dXh1 -= boxX * rintf((dXh1 * boxXinv));
       float dYh1 =  pY - h1Y;
       dYh1 -= boxY * rintf((dYh1 * boxYinv));
       float dZh1 =  pZ - h1Z;
       dZh1 -= boxZ * rintf((dZh1 * boxZinv));
       float rWh1P2 = dXh1*dXh1 + dYh1*dYh1 + dZh1*dZh1;
       float rWh1P_inv = rsqrtf(rWh1P2);
       float enWh1 = h1Q*pQ*rWh1P_inv;

       // h2
       float dXh2 =  pX - h2X;
       dXh2 -= boxX * rintf((dXh2 * boxXinv));
       float dYh2 =  pY - h2Y;
       dYh2 -= boxY * rintf((dYh2 * boxYinv));
       float dZh2 =  pZ - h2Z;
       dZh2 -= boxZ * rintf((dZh2 * boxZinv));
       float rWh2P2 = dXh2*dXh2 + dYh2*dYh2 + dZh2*dZh2;
       float rWh2P_inv = rsqrtf(rWh2P2);
       float enWh2 = h2Q*pQ*rWh2P_inv;

      en_Prot+=(enWo + enWh1 + enWh2);

     } // end for
      dVFull_dev[idx_water] = en_Prot;

      float en_dV = 0.0f;
      for(int i=0;i<NdV;i++)
      {
       float4 dvos_pos = dvos_dev[i];
       float pX = dvos_pos.x;
       float pY = dvos_pos.y;
       float pZ = dvos_pos.z; 
       float pQ = dvos_pos.w;

       // oxygen
       float dXo =  pX - oX;
       dXo -= boxX * rintf((dXo * boxXinv));
       float dYo =  pY - oY;
       dYo -= boxY * rintf((dYo * boxYinv));
       float dZo =  pZ - oZ;
       dZo -= boxZ * rintf((dZo * boxZinv));
       float rWoP2 = dXo*dXo + dYo*dYo + dZo*dZo;
       float rWoP_inv = rsqrtf(rWoP2);
       float enWo = oQ*pQ*rWoP_inv;

       // h1
       float dXh1 =  pX - h1X;
       dXh1 -= boxX * rintf((dXh1 * boxXinv));
       float dYh1 =  pY - h1Y;
       dYh1 -= boxY * rintf((dYh1 * boxYinv));
       float dZh1 =  pZ - h1Z;
       dZh1 -= boxZ * rintf((dZh1 * boxZinv));
       float rWh1P2 = dXh1*dXh1 + dYh1*dYh1 + dZh1*dZh1;
       float rWh1P_inv = rsqrtf(rWh1P2);
       float enWh1 = h1Q*pQ*rWh1P_inv;

       // h2
       float dXh2 =  pX - h2X;
       dXh2 -= boxX * rintf((dXh2 * boxXinv));
       float dYh2 =  pY - h2Y;
       dYh2 -= boxY * rintf((dYh2 * boxYinv));
       float dZh2 =  pZ - h2Z;
       dZh2 -= boxZ * rintf((dZh2 * boxZinv));
       float rWh2P2 = dXh2*dXh2 + dYh2*dYh2 + dZh2*dZh2;
       float rWh2P_inv = rsqrtf(rWh2P2);
       float enWh2 = h2Q*pQ*rWh2P_inv;

       en_dV+=(enWo + enWh1 + enWh2);

     } // end for

       dVDelta_dev[idx_water] = en_dV;

     }
   __syncthreads();
     return; 



}
/*
   cosShell_GPU<<<numBlocks, numThreadsPerBlock>>>(nSh,*prMx,*prMy,*prMz,*dCos
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  water_binfo_dev);


*/
__global__ void cosShell_GPU(int nSh,float Px,float Py,float Pz, float dCos,
                                                  float4* water_oxy_dev,float4* water_h1_dev,
                                                  float4* water_h2_dev,
                                                  float4* water_binfo_dev)
{
     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < nSh)
     {// then we actually do something...
  
      // get our oxygen's atom's position
      float4 water_o_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float oX = water_o_pos.x;
      float oY = water_o_pos.y;
      float oZ = water_o_pos.z;
      float oQ = water_o_pos.w; // fortran's water index (float)

   // get our hydrogen's #1 atom's position
      float4 water_h1_pos = water_h1_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h1X = water_h1_pos.x;
      float h1Y = water_h1_pos.y;
      float h1Z = water_h1_pos.z;
      float h1Q = water_h1_pos.w; // fortran's water index (float)

   // get our hydrogen's #1 atom's position
      float4 water_h2_pos = water_h2_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h2X = water_h2_pos.x;
      float h2Y = water_h2_pos.y;
      float h2Z = water_h2_pos.z;
      float h2Q = water_h2_pos.w; // fortran's water index (float)

  // find the waters dipole moment
      float eX=oX*oQ;
      float eY=oY*oQ;
      float eZ=oZ*oQ;

      eX+=h1X*h1Q;
      eY+=h1Y*h1Q;
      eZ+=h1Z*h1Q;

      eX+=h2X*h2Q;
      eY+=h2Y*h2Q;
      eZ+=h2Z*h2Q;
 
      float e2=eX*eX + eY*eY + eZ*eZ;
      float elen_inv = rsqrtf(e2);
      eX*=elen_inv; 
      eY*=elen_inv; 
      eZ*=elen_inv; 

      float cosEW = Px*eX + Py*eY + Pz*eZ;

      if( cosEW < -1.0f) cosEW = -1.0f;
      if( cosEW > 1.0f) cosEW = 1.0f;
  
      float binCEW = cosEW/dCos;

      int ibinCEW = __float2int_rn(binCEW);
      ibinCEW++;
      water_binfo_dev[idx_water].x = cosEW;
      water_binfo_dev[idx_water].y = __float2int_rn(ibinCEW);

     }
   __syncthreads();
     return; 

}
/*
     
   calcProtMu_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtAtoms,*ddX,*ddY,*ddZ,
                                                     protein_dev,protMu_dev);

*/
__global__ void calcProtMu_GPU(int nP,float mX,float mY,float mZ, float4* protein_dev,float4* protMu_dev)
{
     unsigned int idx_prot = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_prot < nP)
     {// then we actually do something...

  
      // get our oxygen's atom's position
      float4 prot_pos = protein_dev[idx_prot];//tex1Dfetch(water_tex, idx_water);   
      float pX = prot_pos.x;
      float pY = prot_pos.y;
      float pZ = prot_pos.z;
      float pM = prot_pos.w; // the mass of the atom

      protMu_dev[idx_prot].x = (pX - mX)*pM;
      protMu_dev[idx_prot].y = (pY - mY)*pM;
      protMu_dev[idx_prot].z = (pZ - mZ)*pM;
      protMu_dev[idx_prot].w = pM;

     }              
   __syncthreads();
     return; 

}

/*

   p2shell_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NProtAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,protein_dev,
                                                  p2arg_dev);


*/
__global__ void p2shell_GPU(float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int Nw, int Np, bool isTruncated, float4* water_oxy_dev, float4* water_h1_dev, float4* water_h2_dev, float4* protein_dev, float* p2arg_dev)
{

     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < Nw)
     {// then we actually do something...

      float rX,rY,rZ;
  
      // get our oxygen's atom's position
      float4 water_o_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float oX = water_o_pos.x;
      float oY = water_o_pos.y;
      float oZ = water_o_pos.z;
      float oQ = water_o_pos.w; // fortran's water index (float)

   // get our hydrogen's #1 atom's position
      float4 water_h1_pos = water_h1_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h1X = water_h1_pos.x;
      float h1Y = water_h1_pos.y;
      float h1Z = water_h1_pos.z;
      float h1Q = water_h1_pos.w; // fortran's water index (float)

   // get our hydrogen's #1 atom's position
      float4 water_h2_pos = water_h2_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   
      float h2X = water_h2_pos.x;
      float h2Y = water_h2_pos.y;
      float h2Z = water_h2_pos.z;
      float h2Q = water_h2_pos.w; // fortran's water index (float)
   

      float eX=oX*oQ;
      float eY=oY*oQ;
      float eZ=oZ*oQ;

      eX+=h1X*h1Q;
      eY+=h1Y*h1Q;
      eZ+=h1Z*h1Q;

      eX+=h2X*h2Q;
      eY+=h2Y*h2Q;
      eZ+=h2Z*h2Q;
 
      float e2=eX*eX + eY*eY + eZ*eZ;
      float elen_inv = rsqrtf(e2);
      eX*=elen_inv; 
      eY*=elen_inv; 
      eZ*=elen_inv; 

      float low_rWP = 99999.;
      for(int i=0;i<Np;i++)
      {
       float4 prot_pos = protein_dev[i];
       float pX = prot_pos.x;
       float pY = prot_pos.y;
       float pZ = prot_pos.z; 

       float dX =  pX - oX;
       dX -= boxX * rintf((dX * boxXinv));

       float dY =  pY - oY;
       dY -= boxY * rintf((dY * boxYinv));

       float dZ =  pZ - oZ;
       dZ -= boxZ * rintf((dZ * boxZinv));

       float rWP2 = dX*dX + dY*dY + dZ*dZ;
       float rWP_inv = rsqrtf(rWP2);
       float rWP = 1.0f/rWP_inv;
      
       if(rWP < low_rWP)
       {
         low_rWP = rWP;
         rX = dX*rWP_inv;
         rY = dY*rWP_inv;
         rZ = dZ*rWP_inv;       
       } //end if

     } // end for
       float pp2 = eX*rX + eY*rY + eZ*rZ;

       p2arg_dev[idx_water] = pp2;
    
   }// else we just sync our threads and go...

   __syncthreads();

     return; 

}


/*

     // check for shell waters based on all water atoms
     isFirstShellAll_GPU: A GPU kernel for determining which waters are first shell waters, where all atoms lie within
                          a given cutoff distance.

     // check for shell waters based on all water atoms
     isFirstShellAll_GPU<<<numBlocks, numThreadsPerBlock>>>(*rCut,
                                     *boxX,*boxY,*boxZ,
                                      boxXinv,boxYinv,boxZinv,
                                      nW3,NProtAtoms,
                                      isTruncated,
                                      water_oxy_dev,
                                      water_h1_dev,water_21_dev,
                                      protein_dev,lshell_dev);


*/
__global__ void isFirstShellAll_GPU(float rC, float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int Nw, int Np, bool isTruncated, float4* water_oxy_dev, float4* water_h1_dev, float4* water_h2_dev, float4* protein_dev, int* lshell_dev)
{

     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < Nw)
     {// then we actually do something...


   // get our oxygen's atom's position
   float4 oxy_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);      
   float oX = oxy_pos.x;
   float oY = oxy_pos.y;
   float oZ = oxy_pos.z;

   float wInd = oxy_pos.w; // fortran's water index (float)
   int wI = __float2int_rn(wInd); // fortran's water index (int)

   // get our h1 atom's position
   float4 h1_pos = water_h1_dev[idx_water];//tex1Dfetch(water_tex, idx_water);      
   float h1X = h1_pos.x;
   float h1Y = h1_pos.y;
   float h1Z = h1_pos.z;

   // get our h2 atom's position
   float4 h2_pos = water_h2_dev[idx_water];//tex1Dfetch(water_tex, idx_water);      
   float h2X = h2_pos.x;
   float h2Y = h2_pos.y;
   float h2Z = h2_pos.z;

   int idx = 0;

    for(int i=0;i<Np;++i)
    {
     float4 prot_pos = protein_dev[i];
     float pX = prot_pos.x;
     float pY = prot_pos.y;
     float pZ = prot_pos.z;
     float pR = prot_pos.w;

     //oxy atom
     float dX =  pX - oX;
     dX -= boxX * rintf((dX * boxXinv));
     float dY =  pY - oY;
     dY -= boxY * rintf((dY * boxYinv));
     float dZ =  pZ - oZ;
     dZ -= boxZ * rintf((dZ * boxZinv));

     float rOP2 = dX*dX + dY*dY + dZ*dZ;
     float rOP = sqrtf(rOP2);

     //h1 atom
     dX =  pX - h1X;
     dX -= boxX * rintf((dX * boxXinv));
     dY =  pY - h1Y;
     dY -= boxY * rintf((dY * boxYinv));
     dZ =  pZ - h1Z;
     dZ -= boxZ * rintf((dZ * boxZinv));

     float rH1P2 = dX*dX + dY*dY + dZ*dZ;
     float rH1P = sqrtf(rH1P2);

     //h2 atom
     dX =  pX - h2X;
     dX -= boxX * rintf((dX * boxXinv));
     dY =  pY - h2Y;
     dY -= boxY * rintf((dY * boxYinv));
     dZ =  pZ - h2Z;
     dZ -= boxZ * rintf((dZ * boxZinv));

     float rH2P2 = dX*dX + dY*dY + dZ*dZ;
     float rH2P = sqrtf(rH2P2);

     // find the acceptable cutoff for this protein atom
     float myCut = pR + rC;

     if((rOP < myCut) && (rH1P < myCut) && (rH2P < myCut))
     {
       //printf("IN KERNEL ::: FOUND cutoff water, rWP = %f, myCut = %f \n", rWP, myCut);
       idx = wI;
       break;
     }

    }
 
    lshell_dev[idx_water] = idx; // update the lshell array with our index or 0

   }// else we just sync our threads and go...

   __syncthreads();

     return; 

}

/*

isFirstShell_GPU: A GPU kernel for determining which waters are actually first shell atoms


*/
__global__ void isFirstShell_GPU(float rC, float boxX, float boxY, float boxZ, float boxXinv, float boxYinv, float boxZinv, int Nw, int Np, bool isTruncated, float4* water_oxy_dev, float4* protein_dev, int* lshell_dev)
{

     unsigned int idx_water = blockIdx.x*blockDim.x + threadIdx.x;

     if (idx_water < Nw)
     {// then we actually do something...


   // get our water's atom's position
    float4 water_pos = water_oxy_dev[idx_water];//tex1Dfetch(water_tex, idx_water);   

   
   float wX = water_pos.x;
   float wY = water_pos.y;
   float wZ = water_pos.z;
   float wInd = water_pos.w; // fortran's water index (float)
   int wI = __float2int_rn(wInd); // fortran's water index (int)

   int idx = 0;

    for(int i=0;i<Np;++i)
    {
     float4 prot_pos = protein_dev[i];
     float pX = prot_pos.x;
     float pY = prot_pos.y;
     float pZ = prot_pos.z;
     float pR = prot_pos.w;

     float dX =  pX - wX;
     dX -= boxX * rintf((dX * boxXinv));

     float dY =  pY - wY;
     dY -= boxY * rintf((dY * boxYinv));

     float dZ =  pZ - wZ;
     dZ -= boxZ * rintf((dZ * boxZinv));

     float rWP2 = dX*dX + dY*dY + dZ*dZ;
     float rWP = sqrtf(rWP2);

     // find the acceptable cutoff for this protein atom
     float myCut = pR + rC;

     if(rWP < myCut)
     {
       //printf("IN KERNEL ::: FOUND cutoff water, rWP = %f, myCut = %f \n", rWP, myCut);
       idx = wI;
       break;
     }

    }
 
    lshell_dev[idx_water] = idx; // update the lshell array with our index or 0

   }// else we just sync our threads and go...

   __syncthreads();

     return; 

}

extern "C" void alloc_gpu_shell_(int* rank, int* iWatStart, int* iPrEnd,int* iAtoms,int* iDVAtoms, int* iAttempts, int* iSelect, int* iGPUs, bool* doGPUInit)
{

   iWaterIdx = (*iWatStart * 3) - 3; // the starting water atomic index
   iLastWat = (*iAtoms - 1) * 3; // the ending water atomic index
   iProtEnd = (*iPrEnd * 3) - 3; // the ending protein atomic index
   iDVEnd = (*iDVAtoms*3) - 3; // the ending DVAtom index
   NDVAtoms = *iDVAtoms; // the number of DVAtoms
   NwAtoms = *iAtoms - *iWatStart + 1;
   NProtAtoms = *iPrEnd;
   NSelAtoms = *iSelect;
   NMCAttempts = *iAttempts;
   iFortWatIdx = *iWatStart;
   // Find our local GPU rank on the node
   //int rank_gpu = *rank % *iGPUs;
   //cudaSetDevice(rank_gpu); // THIS COULD BE A DANGEROUS PLACE TO MAKE THIS CALL

   //char *methodName;
   //methodName = "alloc_gpu_shell";
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
    checkCUDAError2(cudaOut);
    *doGPUInit = false;
   }

   //cudaSetDevice(*rank); // THIS COULD BE A DANGEROUS PLACE TO MAKE THIS CALL

   //char cudaOut[100];
   // int iChr;
   // iChr = sprintf(cudaOut,"cudaSetDevice -- setting device to %i on rank %i \n",rank_gpu,*rank);
   // checkCUDAError2(cudaOut);

   //------------------------------------------
   // PUT THIS IN AN ALLOCATION ROUTINE
   //------------------------------------------
   memWaterF4 = sizeof(float4) *  (NwAtoms/3); 
   memProtF4 = sizeof(float4) *  NProtAtoms;
   memSelF4  = sizeof(float4) *  NSelAtoms;
   memDVAtomsF4 = sizeof(float4) *  NDVAtoms;
   memMCVolF4  = sizeof(float4) * NMCAttempts;
   memMCVolF1  = sizeof(float) * NMCAttempts;

   memShellInt = sizeof(int) *  (NwAtoms/3);
   memFloatWat = sizeof(float) * (NwAtoms/3);

   mc_ran_host = (float4*)malloc(memMCVolF4);
   mc_prot_host = (float*)malloc(memMCVolF1);
   mc_wat_host = (float*)malloc(memMCVolF1);
   mc_count_host = (float*)malloc(memMCVolF1);

   water_oxy_host = (float4*)malloc(memWaterF4);
   water_h1_host  = (float4*)malloc(memWaterF4);
   water_h2_host  = (float4*)malloc(memWaterF4);
   water_binfo_host = (float4*)malloc(memWaterF4);

   protein_host   = (float4*)malloc(memProtF4);
   protein_sel_host   = (float4*)malloc(memSelF4);
   protMu_host   = (float4*)malloc(memProtF4);
   protMu_sel_host   = (float4*)malloc(memSelF4);

   dvos_host  = (float4*)malloc(memDVAtomsF4);

   p2arg_host     = (float*)malloc(memFloatWat);
   lshell_host    = (int*)malloc(memShellInt);

   dVpol_host  = (float*)malloc(memFloatWat);

   dQpol_host  = (float*)malloc(memFloatWat);

   dVFull_host = (float*)malloc(memFloatWat);
   dVDelta_host = (float*)malloc(memFloatWat);

   water_ef_host = (float*)malloc(memFloatWat);

   cudaMalloc( (void **) &water_oxy_dev, memWaterF4);
  // checkCUDAError2("cudaMalloc -- water_oxy_dev");
   cudaMalloc( (void **) &water_h1_dev, memWaterF4);
  // checkCUDAError2("cudaMalloc -- water_h1_dev");
   cudaMalloc( (void **) &water_h2_dev, memWaterF4);
   cudaMalloc( (void **) &water_binfo_dev, memWaterF4);
  // checkCUDAError2("cudaMalloc -- water_h2_dev");

   cudaMalloc( (void **) &protein_dev,memProtF4);
   cudaMalloc( (void **) &protein_sel_dev,memSelF4);
  // checkCUDAError2("cudaMalloc -- protein_dev");
   cudaMalloc( (void **) &protMu_dev,memProtF4);
   cudaMalloc( (void **) &protMu_sel_dev,memSelF4);
  // checkCUDAError2("cudaMalloc -- protMu_dev");

   cudaMalloc( (void **) &dvos_dev,memDVAtomsF4);
  // checkCUDAError2("cudaMalloc -- dvos_dev");

   cudaMalloc( (void **) &mc_ran_dev,memMCVolF4);
  // checkCUDAError2("cudaMalloc -- mc_ran_dev");
   cudaMalloc( (void **) &mc_prot_dev,memMCVolF1);
  // checkCUDAError2("cudaMalloc -- mc_ran_dev");
   cudaMalloc( (void **) &mc_wat_dev,memMCVolF1);
  // checkCUDAError2("cudaMalloc -- mc_wat_dev");
   cudaMalloc( (void **) &mc_count_dev,memMCVolF1);
  // checkCUDAError2("cudaMalloc -- mc_count_dev");

   cudaMalloc( (void **) &lshell_dev,memShellInt);
  // checkCUDAError2("cudaMalloc -- lshell_dev");

   cudaMalloc( (void **) &p2arg_dev, memFloatWat);
  // checkCUDAError2("cudaMalloc -- p2arg_dev");

   cudaMalloc( (void **) &dVpol_dev, memFloatWat);
  // checkCUDAError2("cudaMalloc -- dVpol_dev");

   cudaMalloc( (void **) &dQpol_dev, memFloatWat);
  // checkCUDAError2("cudaMalloc -- dQpol_dev");

   cudaMalloc( (void **) &dVFull_dev, memFloatWat);
  // checkCUDAError2("cudaMalloc -- dVFull_dev");
   cudaMalloc( (void **) &dVDelta_dev, memFloatWat);
  // checkCUDAError2("cudaMalloc -- dVDelta_dev");

   cudaMalloc( (void **) &water_ef_dev,memFloatWat);
  // checkCUDAError2("cudaMalloc -- water_ef_dev");

   //printf("ALLOCATED ALL HOST/GPU memory\n");
   //------------------------------------------
   // END ALLOCATION ROUTINE
   //------------------------------------------

}

extern "C" void dealloc_gpu_shell_()
{
   //------------------------------------------
   // PUT THIS IN A DEALLOCATION ROUTINE
   //------------------------------------------
   // free device memory

    cudaFree(water_ef_dev);
    cudaFree(water_oxy_dev);
    cudaFree(water_h1_dev);
    cudaFree(water_h2_dev);
    cudaFree(water_binfo_dev);

    cudaFree(protein_dev);
    cudaFree(protein_sel_dev);
    cudaFree(protMu_dev);
    cudaFree(protMu_sel_dev);

    cudaFree(dvos_dev);

    cudaFree(lshell_dev);
    cudaFree(p2arg_dev);

    cudaFree(dVFull_dev);
    cudaFree(dVDelta_dev);

    cudaFree(dVpol_dev);
    cudaFree(dQpol_dev);

    cudaFree(mc_ran_dev);
    cudaFree(mc_prot_dev);
    cudaFree(mc_wat_dev);
    cudaFree(mc_count_dev);

   // free up host memory
      delete[] dVpol_host;
      dVpol_host = NULL;

      delete[] dQpol_host;
      dQpol_host = NULL;

      delete[] dVFull_host;
      dVFull_host = NULL;
      delete[] dVDelta_host;
      dVDelta_host = NULL;

      delete[] p2arg_host;
      p2arg_host = NULL;

      delete[] dvos_host;
      dvos_host = NULL;

      delete[] water_binfo_host;
      water_binfo_host = NULL;

      delete[] water_oxy_host;
      water_oxy_host = NULL;

      delete[] water_h1_host;
      water_h1_host = NULL;

      delete[] water_h2_host;
      water_h2_host = NULL;
     
      delete[] protein_host;
      protein_host = NULL;

      delete[] protein_sel_host;
      protein_sel_host = NULL;

      delete[] protMu_host;
      protMu_host = NULL;

      delete[] protMu_sel_host;
      protMu_sel_host = NULL;

      delete[] lshell_host;
      lshell_host = NULL;

      delete[] water_ef_host;
      water_ef_host = NULL;

      delete[] mc_ran_host;
      mc_ran_host = NULL;

      delete[] mc_prot_host;
      mc_prot_dev = NULL;

      delete[] mc_wat_host;
      mc_wat_host = NULL;

      delete[] mc_count_host;
      mc_count_host = NULL;
   //------------------------------------------
   // END DEALLOCATION ROUTINE
   //------------------------------------------

}

/*

get_shell_gpu_: (called from fortran as "get_shell_gpu") is the fortran interface function for the isFirstShell_GPU
                kernel code


*/
extern "C" void get_shell_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* rCut, int* lShell, int* nShell, bool* checkForAll)
{
   bool isTruncated = false;
   int nW3 = NwAtoms/3;

   // setup the grid to run the kernel
   int numThreadsPerBlock = 256;
   //int numThreadsPerBlock = 1;//256/8;
   int numBlocks = nW3 / numThreadsPerBlock + (nW3 % numThreadsPerBlock == 0?0:1);

   float* xx;
   float* yy;
   float* zz;

   float* h1x;
   float* h1y;
   float* h1z;

   float* h2x;
   float* h2y;
   float* h2z;

   //float* rr;
   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);


   // Start at the first index for the protein
   xx = &xyz[0];
   yy = &xyz[1];
   zz = &xyz[2];
   //rr = &qrms[1]; /// points to 1 for the radius position
   int j = -3;

   // loop over all protein atoms and fill up the protein array
   for(int i = 0;i < NProtAtoms;i++) 
   {
    j+=3;
    // fill up the protein on the host in float4s
    //protein_host[i] = make_float4(xx[j],yy[j],zz[j],rr[j]);
    protein_host[i] = make_float4(xx[j],yy[j],zz[j],0.0f);
   }

   xx = &xyz[iWaterIdx];
   yy = &xyz[iWaterIdx+1];
   zz = &xyz[iWaterIdx+2];
   j = -9; // 9 = 3atoms*(1x+1y+1z)
   float watNum;
   // loop over all oxygens in the waters and fill up the water array
   for(int ii = 0;ii < nW3;ii++)
   {
    j+=9;
    watNum = float(iFortWatIdx + 3*ii);
    water_oxy_host[ii] = make_float4(xx[j],yy[j],zz[j],watNum);
   }

    // copy the water from the host to the device
   cudaMemcpy(protein_dev, protein_host, memProtF4, cudaMemcpyHostToDevice);
   checkCUDAError2("cudaMemcpy(prot) -- host to device in function get_shell_gpu");
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
   checkCUDAError2("cudaMemcpy(wat) -- host to device in function get_shell_gpu");

   if(*checkForAll) 
   {

     // pointers to the h1 atom
     h1x = &xyz[iWaterIdx+3];
     h1y = &xyz[iWaterIdx+4];
     h1z = &xyz[iWaterIdx+5];

     // pointers to the h2 atom
     h2x = &xyz[iWaterIdx+6];
     h2y = &xyz[iWaterIdx+7];
     h2z = &xyz[iWaterIdx+8];

     j = -9; // 9 = 3atoms*(1x+1y+1z)
     float watNum;
     // loop over all oxygens in the waters and fill up the water array
     for(int ii = 0;ii < nW3;ii++)
     {
      j+=9;
      watNum = float(iFortWatIdx + 3*ii);
      water_h1_host[ii] = make_float4(h1x[j],h1y[j],h1z[j],watNum);
      water_h2_host[ii] = make_float4(h2x[j],h2y[j],h2z[j],watNum);
     }

     cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
     cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);


     // check for shell waters based on all water atoms
     isFirstShellAll_GPU<<<numBlocks, numThreadsPerBlock>>>(*rCut,
                                     *boxX,*boxY,*boxZ,
                                      boxXinv,boxYinv,boxZinv,
                                      nW3,NProtAtoms,
                                      isTruncated,
                                      water_oxy_dev,
                                      water_h1_dev,water_h2_dev,
                                      protein_dev,lshell_dev);

   } else {
     // only check for shell waters based on the O-atom only
     isFirstShell_GPU<<<numBlocks, numThreadsPerBlock>>>(*rCut,
                                     *boxX,*boxY,*boxZ,
                                      boxXinv,boxYinv,boxZinv,
                                      nW3,NProtAtoms,
                                         isTruncated,
                         water_oxy_dev,protein_dev,lshell_dev);


   }
//checkCUDAError2("isFirstShell_GPU -- KERNEL ERROR!!!");
  // copy the water from the host to the device
   cudaMemcpy(lshell_host,lshell_dev,memShellInt,cudaMemcpyDeviceToHost);
  // checkCUDAError2("cudaMemcpy -- results from device to host");


    
  // Sum the whole array lshell_host array
    int myN = 0;
    for(int i=0;i<nW3;i++)
    {
      
      if(lshell_host[i] > 0)
      {
      // if(myN==100)printf("IN CUDA::: shell val[%i] = %i \n",myN,lshell_host[i]);
       lShell[myN] = lshell_host[i];
        myN++;
      }

    }
    *nShell = myN;
    //printf(" IN CUDA::: myN = %i \n",myN);
}


/*

get_p2_gpu_: (called from fortran as "get_p2_gpu") is the fortran interface function for the 
                   shellp2 kernel on the GPU.  
                   NOTE: THIS KERNEL ASSUMES THAT PROTEIN HAS BEEN COPIED TO THE GPU

*/
extern "C" void get_p2_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, int* lShell,int* nShell,float* p2arg, float* p2Shell)
{
   bool isTruncated = false;
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = nSh / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float* xh1;
   float* yh1;
   float* zh1;

   float* xh2;
   float* yh2;
   float* zh2;

   float* qo;
   float* qh1;
   float* qh2;

   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);


   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];
   qo = &qrms[iWatShell]; // oxygen

   xh1 = &xyz[iWatShell+3];
   yh1 = &xyz[iWatShell+4];
   zh1 = &xyz[iWatShell+5];
   qh1 = &qrms[iWatShell+3]; // h1

   xh2 = &xyz[iWatShell+6];
   yh2 = &xyz[iWatShell+7];
   zh2 = &xyz[iWatShell+8];
   qh2 = &qrms[iWatShell+6]; // h2

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,*qo);
    water_h1_host[i] = make_float4(*xh1,*yh1,*zh1,*qh1);
    water_h2_host[i] = make_float4(*xh2,*yh2,*zh2,*qh2);

   }

     // copy the water from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
   //checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
   //checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);
   //checkCUDAError2("cudaMemcpy(wat) -- host to device");

   p2shell_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NProtAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,protein_dev,
                                                  p2arg_dev);

   // copy the p2args from the host to the device
   cudaMemcpy(p2arg_host,p2arg_dev,memFloatWat,cudaMemcpyDeviceToHost);
   //checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum the whole array p2arg_host array, and calculate the order parameter
    float p2 = 0.0f;
    float p2Sh = 0.0f;
    float parg;
    for(int i=0;i<nSh;i++)
    {
      parg=p2arg_host[i];
      p2+=parg;
      p2Sh+=(3.0f*(parg*parg) - 1.0f);
    }
    *p2arg = p2;
    *p2Shell = p2Sh;





}
/*

get_dvos_shl_gpu_: (called from fortran as "get_dv0s_shl_gpu") is the fortran interface function for the 
                   dvos_shl_GPU kernel on the GPU.  
                   NOTE: THIS KERNEL ASSUMES THAT THE PROTEIN HAS BEEN COPIED OVER

*/
extern "C" void get_dv0s_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvAtoms, float* dvQs, int* lShell,int* nShell,float* dVFull, float* dVDelta)
{
   bool isTruncated = false;
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = nSh / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float* xh1;
   float* yh1;
   float* zh1;

   float* xh2;
   float* yh2;
   float* zh2;

   float* qo;
   float* qh1;
   float* qh2;

   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);


   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];
   qo = &qrms[iWatShell]; // oxygen

   xh1 = &xyz[iWatShell+3];
   yh1 = &xyz[iWatShell+4];
   zh1 = &xyz[iWatShell+5];
   qh1 = &qrms[iWatShell+3]; // h1

   xh2 = &xyz[iWatShell+6];
   yh2 = &xyz[iWatShell+7];
   zh2 = &xyz[iWatShell+8];
   qh2 = &qrms[iWatShell+6]; // h2

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,*qo);
    water_h1_host[i] = make_float4(*xh1,*yh1,*zh1,*qh1);
    water_h2_host[i] = make_float4(*xh2,*yh2,*zh2,*qh2);

   }

   float* qq = &qrms[0]; /// points to 1 for the radius position
   int j = -3;
   // loop over all protein atoms and fill up the protein array
   for(int i = 0;i < NProtAtoms;i++) 
   {
    j+=3;
    // fill up the protein charges the host in the w position of the float4s
    protein_host[i].w = qq[j];
   }

   float* dq = &dvQs[0];
   float* xdq = &dvAtoms[0];
   float* ydq = &dvAtoms[1];
   float* zdq = &dvAtoms[2];

//   float tmpQ = 0.0f;
   j = -3;
   // loop over all dv0s atoms and fill up the dv0s array
   for(int i = 0;i < NDVAtoms;i++) 
   {
    j+=3;
    // fill up the protein charges the host in the w position of the float4s
//    tmpQ+=dq[i];
    dvos_host[i].x = xdq[j]; 
    dvos_host[i].y = ydq[j]; 
    dvos_host[i].z = zdq[j]; 
    dvos_host[i].w = dq[i];
   }

//    printf("IN CPU::: total dQ = %f for %i atoms \n",tmpQ,NDVAtoms);

    // copy the dv0s atoms (with charges) from the host to the device
   cudaMemcpy(dvos_dev, dvos_host, memDVAtomsF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(dvos) -- host to device");
    // copy the protein (with charges) from the host to the device
   cudaMemcpy(protein_dev, protein_host, memProtF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(protein) -- host to device");
     // copy the water(with charges) from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");

   dvos_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NProtAtoms,NDVAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  protein_dev,dvos_dev,
                                                  dVFull_dev,dVDelta_dev);
//   checkCUDAError2("dvos_shl_GPU");

   // copy the dV0s arrays from the host to the device
   cudaMemcpy(dVFull_host,dVFull_dev,memFloatWat,cudaMemcpyDeviceToHost);
   cudaMemcpy(dVDelta_host,dVDelta_dev,memFloatWat,cudaMemcpyDeviceToHost);

//   checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum both dVFull/dVDelta arrays, and calculate the full energies
    float DV = 0.0f;
    float DVDel = 0.0f;
    for(int i=0;i<nSh;i++)
    {
      DV+=dVFull_host[i];
      DVDel+=dVDelta_host[i];

//      printf("IN CPU:: dVFull_host[%i] = %f   dVDelta_host[i] = %f \n",i,dVFull_host[i],i,dVDelta_host[i]);
    }

//    printf("IN CPU:: DV = %f   DVDel = %f \n",DV,DVDel);
    *dVFull = DV;
    *dVDelta = DVDel;

}

/*

get_vol_shl_gpu_: (called from fortran as "get_vol_shl_gpu") is the fortran interface function for the 
                   vol_shl_GPU kernel on the GPU.  
                   NOTE: THIS KERNEL ASSUMES THAT THE PROTEIN HAS BEEN COPIED OVER

*/
extern "C" void get_vol_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, int* lShell,int* nShell, float* mcRadius, float* randx, float* randy, float* randz, float* centerX, float* centerY, float* centerZ, float* volProt, float* volWat, float* epsProt, float* epsWat)
{
   bool isTruncated = false;
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = NMCAttempts / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float* xh1;
   float* yh1;
   float* zh1;

   float* xh2;
   float* yh2;
   float* zh2;

   float* ro;
   float* rh1;
   float* rh2;

   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);

   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];
   ro = &qrms[iWatShell+1]; // oxygen

   xh1 = &xyz[iWatShell+3];
   yh1 = &xyz[iWatShell+4];
   zh1 = &xyz[iWatShell+5];
   rh1 = &qrms[iWatShell+4]; // h1



   xh2 = &xyz[iWatShell+6];
   yh2 = &xyz[iWatShell+7];
   zh2 = &xyz[iWatShell+8];
   rh2 = &qrms[iWatShell+7]; // h2

//   if(i % 100 == 0) printf("h2 radius = %f \n",rh2);

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,*ro);
    water_h1_host[i] = make_float4(*xh1,*yh1,*zh1,*rh1);
    water_h2_host[i] = make_float4(*xh2,*yh2,*zh2,*rh2);

   }

   float* rr = &qrms[1]; /// points to 1 for the radius position
   int j = -3;
   // loop over all protein atoms and fill up the protein array
   for(int i = 0;i < NProtAtoms;i++) 
   {
    j+=3;
    // fill up the protein charges the host in the w position of the float4s
    protein_host[i].w = rr[j];
   }

   for(int i = 0;i < NMCAttempts;i++)
   {
    mc_ran_host[i].x = randx[i];
    mc_ran_host[i].y = randy[i];
    mc_ran_host[i].z = randz[i];

//   if(i % 100 == 0) printf("rx,ry,rz = %f %f %f \n",randx[i],randy[i],randz[i]);
   }



   // copy all random numbers to the GPU
   //    ---> seriously man, move this to the GPU
   cudaMemcpy(mc_ran_dev, mc_ran_host, memMCVolF4, cudaMemcpyHostToDevice);


    // copy the protein (with radii) from the host to the device
   cudaMemcpy(protein_dev, protein_host, memProtF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(protein) -- host to device");
     // copy the water(with charges) from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");

   vol_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  *mcRadius,
                                                  *centerX, *centerY, *centerZ,
                                                  nSh,NProtAtoms,NMCAttempts,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,water_h2_dev,
                                                  protein_dev,
                                                  mc_ran_dev,
                                                  mc_prot_dev,mc_wat_dev,mc_count_dev);
//   checkCUDAError2("dvos_shl_GPU");

   // copy the dV0s arrays from the host to the device
   cudaMemcpy(mc_prot_host,mc_prot_dev,memMCVolF1,cudaMemcpyDeviceToHost);
   cudaMemcpy(mc_wat_host,mc_wat_dev,memMCVolF1,cudaMemcpyDeviceToHost);
   cudaMemcpy(mc_count_host,mc_count_dev,memMCVolF1,cudaMemcpyDeviceToHost);

//   checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum both dVFull/dVDelta arrays, and calculate the full energies
    float sum_prot = 0.0f;
    float sum_wat = 0.0f;
    float sum_tot = 0.0f;

    float prot_e = 1.0f;
    float prot_ee;
    float wat_e = 1.0f;
    float wat_ee;    

    float epsP, epsW;
    for(int i=0;i<NMCAttempts;i++)
    {
      sum_tot +=  mc_count_host[i];
      sum_prot += mc_prot_host[i];
      sum_wat += mc_wat_host[i];

      if ( i % 10000 == 0)
      {

      //  printf("mc_count_host[i] = %f \n",mc_count_host[i]);
      //  printf("mc_prot_host[i] = %f \n",mc_prot_host[i]);
      //  printf("mc_wat_host[i] = %f \n",mc_wat_host[i]);

        prot_ee = sum_prot/sum_tot;
        epsP = fabs( (prot_e - prot_ee)/prot_ee );
        prot_e  = prot_ee;

        wat_ee = sum_wat/sum_tot;
        epsW = fabs( (wat_e - wat_ee)/wat_ee );
        wat_e  = wat_ee;
      }
//      printf("IN CPU:: dVFull_host[%i] = %f   dVDelta_host[i] = %f \n",i,dVFull_host[i],i,dVDelta_host[i]);
    }

     float sphereVol = (4.0f/3.0f)*3.14159f;
     sphereVol *= *mcRadius;
     sphereVol *= *mcRadius;
     sphereVol *= *mcRadius;

     float vol_prot = (sum_prot/sum_tot)*sphereVol;
     float vol_wat = (sum_wat/sum_tot)*sphereVol;

//    printf("IN CPU:: DV = %f   DVDel = %f \n",DV,DVDel);
    *volProt = vol_prot;
    *volWat = vol_wat;
    *epsProt = epsP;
    *epsWat = epsW;
}

/*

get_vpol_shl_gpu_: (called from fortran as "get_vpol_shl_gpu") is the fortran interface function for the 
                   vpol_shl_GPU kernel on the GPU.  
                   NOTE: THIS KERNEL ASSUMES THAT THE PROTEIN HAS BEEN COPIED OVER

*/
extern "C" void get_vpol_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvAtoms, float* dvQs, int* lShell,int* nShell,float* dVpol)
{
   bool isTruncated = false;
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = nSh / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float* xh1;
   float* yh1;
   float* zh1;

   float* xh2;
   float* yh2;
   float* zh2;

   float* qo;
   float* qh1;
   float* qh2;

   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);


   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];
   qo = &qrms[iWatShell]; // oxygen

   xh1 = &xyz[iWatShell+3];
   yh1 = &xyz[iWatShell+4];
   zh1 = &xyz[iWatShell+5];
   qh1 = &qrms[iWatShell+3]; // h1

   xh2 = &xyz[iWatShell+6];
   yh2 = &xyz[iWatShell+7];
   zh2 = &xyz[iWatShell+8];
   qh2 = &qrms[iWatShell+6]; // h2

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,*qo);
    water_h1_host[i] = make_float4(*xh1,*yh1,*zh1,*qh1);
    water_h2_host[i] = make_float4(*xh2,*yh2,*zh2,*qh2);

   }

   float* dq = &dvQs[0];
   float* xdq = &dvAtoms[0];
   float* ydq = &dvAtoms[1];
   float* zdq = &dvAtoms[2];

//   float tmpQ = 0.0f;
   int j = -3;
   // loop over all dv0s atoms and fill up the dv0s array
   for(int i = 0;i < NDVAtoms;i++) 
   {
    j+=3;
    // fill up the protein charges the host in the w position of the float4s
//    tmpQ+=dq[i];
    dvos_host[i].x = xdq[j]; 
    dvos_host[i].y = ydq[j]; 
    dvos_host[i].z = zdq[j]; 
    dvos_host[i].w = dq[i];
   }

//    printf("IN CPU::: total dQ = %f for %i atoms \n",tmpQ,NDVAtoms);

    // copy the dv0s atoms (with charges) from the host to the device
   cudaMemcpy(dvos_dev, dvos_host, memDVAtomsF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(dvos) -- host to device");
     // copy the water(with charges) from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");

   vpol_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NDVAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  dvos_dev,
                                                  dVpol_dev);
//   checkCUDAError2("dvos_shl_GPU");

   // copy the dVpol array from the host to the device
   cudaMemcpy(dVpol_host,dVpol_dev,memFloatWat,cudaMemcpyDeviceToHost);

//   checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum the DVpol array, and calculate the full Vpol correction
    float DVpol = 0.0f;
    for(int i=0;i<nSh;i++)
    {
      DVpol+=dVpol_host[i];

//      printf("IN CPU:: dVFull_host[%i] = %f   dVDelta_host[i] = %f \n",i,dVFull_host[i],i,dVDelta_host[i]);
    }

//    printf("IN CPU:: DV = %f   DVDel = %f \n",DV,DVDel);
    *dVpol = DVpol;

}

/*

get_qpol_shl_gpu_: (called from fortran as "get_qpol_shl_gpu") is the fortran interface function for the 
                   qpol_shl_GPU kernel on the GPU.  
                   NOTE: THIS KERNEL ASSUMES THAT THE PROTEIN HAS BEEN COPIED OVER

*/
extern "C" void get_qpol_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, int* lShell,int* nShell,float* dQpol)
{
   bool isTruncated = false;
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = nSh / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float* xh1;
   float* yh1;
   float* zh1;

   float* xh2;
   float* yh2;
   float* zh2;

   float* qo;
   float* qh1;
   float* qh2;

   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);


   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];
   qo = &qrms[iWatShell]; // oxygen

   xh1 = &xyz[iWatShell+3];
   yh1 = &xyz[iWatShell+4];
   zh1 = &xyz[iWatShell+5];
   qh1 = &qrms[iWatShell+3]; // h1

   xh2 = &xyz[iWatShell+6];
   yh2 = &xyz[iWatShell+7];
   zh2 = &xyz[iWatShell+8];
   qh2 = &qrms[iWatShell+6]; // h2

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,*qo);
    water_h1_host[i] = make_float4(*xh1,*yh1,*zh1,*qh1);
    water_h2_host[i] = make_float4(*xh2,*yh2,*zh2,*qh2);

   }

    // copy the water(with charges) from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");

   qpol_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,
                                                  isTruncated,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  dQpol_dev);
//   checkCUDAError2("dvos_shl_GPU");

   // copy the dVpol array from the host to the device
   cudaMemcpy(dQpol_host,dQpol_dev,memFloatWat,cudaMemcpyDeviceToHost);

//   checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum the DVpol array, and calculate the full Vpol correction
    float DQpol = 0.0f;
    for(int i=0;i<nSh;i++)
    {
      DQpol+=dQpol_host[i];
    }

    *dQpol = DQpol;

}


/*

get_efld_shl_gpu_: (called from fortran as "get_efld_shl_gpu") is the fortran interface function for the 
                   efield_shl_GPU kernel on the GPU.  

*/
extern "C" void get_efld_shl_gpu_(float* boxX, float* boxY, float* boxZ, float *xyz, float* qrms, float* dvAtoms, float* dvQs, int* lShell,int* nShell,float* efield)
{
   bool isTruncated = false;
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = nSh / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float boxXinv, boxYinv, boxZinv;
 
   boxXinv = 1.0f/(*boxX);
   boxYinv = 1.0f/(*boxY);
   boxZinv = 1.0f/(*boxZ);


   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,0.0f);

   }

   float* dq = &dvQs[0];
   float* xdq = &dvAtoms[0];
   float* ydq = &dvAtoms[1];
   float* zdq = &dvAtoms[2];

   int j = -3;
   // loop over all dv0s atoms and fill up the dv0s array
   for(int i = 0;i < NDVAtoms;i++) 
   {
    j+=3;
    // fill up the dvos charges the host in the w position of the float4s
    dvos_host[i].x = xdq[j]; 
    dvos_host[i].y = ydq[j]; 
    dvos_host[i].z = zdq[j]; 
    dvos_host[i].w = dq[i];
   }

//    printf("IN CPU::: total dQ = %f for %i atoms \n",tmpQ,NDVAtoms);

    // copy the dv0s atoms (with charges) from the host to the device
   cudaMemcpy(dvos_dev, dvos_host, memDVAtomsF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(dvos) -- host to device");
    // copy the protein (with charges) from the host to the device
     // copy the water(with charges) from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
//   checkCUDAError2("cudaMemcpy(wat) -- host to device");

   efield_shl_GPU<<<numBlocks, numThreadsPerBlock>>>(*boxX,*boxY,*boxZ,
                                                  boxXinv,boxYinv,boxZinv,
                                                  nSh,NDVAtoms,
                                                  isTruncated,
                                                  water_oxy_dev,dvos_dev,
                                                  water_ef_dev);

//   checkCUDAError2("dvos_shl_GPU");

   // copy the efield array from the host to the device
   cudaMemcpy(water_ef_host,water_ef_dev,memFloatWat,cudaMemcpyDeviceToHost);
//   checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum both dVFull/dVDelta arrays, and calculate the full energies
    float efMag = 0.0f;
    for(int i=0;i<nSh;i++)
    {
      efMag += water_ef_host[i];

    }

//    printf("IN CPU:: efield = %f  \n",efield);
    *efield = efMag/((float)nSh);

}

/*

get_prmu_gpu_: (called from fortran as "get_prmu_gpu") is the fortran interface function for the 
                kernel to calculate the protein's electric or mass dipole moment on the GPU

*/
extern "C" void get_prmu_gpu_(float* qrms, float* ddX,float* ddY, float* ddZ, 
                              float* prX,float* prY, float* prZ, float* totMorQ, int* massOrQ, bool* doOutPos, float* outputPos)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = NProtAtoms / numThreadsPerBlock + (NProtAtoms % numThreadsPerBlock == 0?0:1);

//   bool doMass;

//   if(*massOrQ == 0)
//   {
//     doMass = true;
//   }else
//   {
//     doMass = false;
//   }
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);
   
   float* mOrQ;

   if(*massOrQ == 0)
   {
     mOrQ = &qrms[2]; /// points to 2 for the mass position
   }else if (*massOrQ == 1)
   {
     mOrQ = &qrms[0]; /// points to 0 for the charge position
   }

   //float* qq = &qrms[0];  /// points to 0 for the q position
   //float* mm = &qrms[2]; /// points to 2 for the mass position
   int j = -3;
   // loop over all protein atoms and fill up the protein array

    if(*massOrQ == 0 || *massOrQ == 1) 
    {
      for(int i = 0;i < NProtAtoms;i++) 
      {
      j+=3;
      // fill up the proteins q/m the host in the w position of the float4s
      protein_host[i].w = mOrQ[j];
      }

    } else {
      for(int i = 0;i < NProtAtoms;i++) 
      {
       protein_host[i].w = 1.000000000;
      }
    }

   
    // copy the water from the host to the device
   cudaMemcpy(protein_dev, protein_host, memProtF4, cudaMemcpyHostToDevice);

   calcProtMu_GPU<<<numBlocks, numThreadsPerBlock>>>(NProtAtoms,*ddX,*ddY,*ddZ,
                                                     protein_dev,protMu_dev);


   // copy the dipoles from the host to the device
   cudaMemcpy(protMu_host,protMu_dev,memProtF4,cudaMemcpyDeviceToHost);
   //checkCUDAError2("cudaMemcpy -- results from device to host");
    
   // Sum the whole array protMu_host array, and calculate the output vector
    float prMX = 0.0f;
    float prMY = 0.0f;
    float prMZ = 0.0f;
    float prMorQ = 0.0f;

    for(int i=0;i<NProtAtoms;i++)
    {
       prMX += protMu_host[i].x;
       prMY += protMu_host[i].y;
       prMZ += protMu_host[i].z;
       prMorQ += protMu_host[i].w;
    }

    if(*doOutPos) 
    {

      float* xx;
      float* yy;
      float* zz;
      xx = &outputPos[0];
      yy = &outputPos[1];
      zz = &outputPos[2];

      j = -3;
      for(int i=0;i<NProtAtoms;i++)
      {
       j+=3;
       xx[j] = protMu_host[i].x;
       yy[j] = protMu_host[i].y;
       zz[j] = protMu_host[i].z;
       }
    }

    *prX = prMX;
    *prY = prMY;
    *prZ = prMZ;
    *totMorQ = prMorQ;
}

/*

get_prmu_sel_gpu_: (called from fortran as "get_prmu_sel_gpu") is the fortran interface function for the 
                kernel to calculate the protein's electric or mass dipole moment (for selected atoms only) on the GPU

*/
extern "C" void get_prmu_sel_gpu_(float* qrms, float* ddX,float* ddY, float* ddZ, 
                              float* prX,float* prY, float* prZ, float* totMorQ, int* massOrQ, int* selIndices, bool* doOutPos, float* outputPos)
{

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 32; // this is small because of the potentially small number of sel atoms
   int numBlocks = NSelAtoms / numThreadsPerBlock + (NSelAtoms % numThreadsPerBlock == 0?0:1);

/* N              = 1,2,3,4, 5, 6, 7, 8,9,10,77*/
/* WANTED         = 0,3,6,9,12,15,18,21,..3N-1 */
/* 3*(N  + 1) - 6 = 

                    0,3,6,9,12 */
   
   float* mOrQ;

   if(*massOrQ == 0)
   {
     mOrQ = &qrms[2]; /// points to 2 for the mass position
   }else if (*massOrQ == 1)
   {
     mOrQ = &qrms[0]; /// points to 2 for the charge position
   }

   //float* qq = &qrms[0];  /// points to 0 for the q position
   //float* mm = &qrms[2]; /// points to 2 for the mass position
   int j,k;
   // loop over all protein atoms and fill up the protein array

    if(*massOrQ == 0 || *massOrQ == 1) 
    {
      for(int i = 0;i < NSelAtoms;i++) 
      {
       k = selIndices[i];
       j = 3*k - 3;

      // fill up the proteins q/m the host in the w position of the float4s
      protein_sel_host[i].w = mOrQ[j];
      }

    } else {
      for(int i = 0;i < NSelAtoms;i++) 
      {
       protein_sel_host[i].w = 1.000000000;
      }
    }

   
    // copy the water from the host to the device
   cudaMemcpy(protein_sel_dev, protein_sel_host, memSelF4, cudaMemcpyHostToDevice);

   calcProtMu_GPU<<<numBlocks, numThreadsPerBlock>>>(NSelAtoms,*ddX,*ddY,*ddZ,
                                                     protein_sel_dev,protMu_sel_dev);


   // copy the dipoles from the host to the device
   cudaMemcpy(protMu_sel_host,protMu_sel_dev,memSelF4,cudaMemcpyDeviceToHost);
   //checkCUDAError2("cudaMemcpy -- results from device to host");
    
   // Sum the whole array p2arg_host array, and calculate the order parameter
    float prMX = 0.0f;
    float prMY = 0.0f;
    float prMZ = 0.0f;
    float prMorQ = 0.0f;

    for(int i=0;i<NSelAtoms;i++)
    {
       prMX += protMu_sel_host[i].x;
       prMY += protMu_sel_host[i].y;
       prMZ += protMu_sel_host[i].z;
       prMorQ += protMu_sel_host[i].w;
    }

    if(*doOutPos) 
    {

      float* xx;
      float* yy;
      float* zz;
      xx = &outputPos[0];
      yy = &outputPos[1];
      zz = &outputPos[2];

      j = -3;
      for(int i=0;i<NSelAtoms;i++)
      {
       j+=3;
       xx[j] = protMu_sel_host[i].x;
       yy[j] = protMu_sel_host[i].y;
       zz[j] = protMu_sel_host[i].z;
       }
    }

    *prX = prMX;
    *prY = prMY;
    *prZ = prMZ;
    *totMorQ = prMorQ;
}

/*

set_prot_gpu_: (called from fortran as "set_prot_gpu") is the fortran interface function for the 
                kernel to set an arbitrary protein structure to protein_dev on the GPU

*/
extern "C" void set_prot_gpu_(float *xyz)
{
    float* xx;
    float* yy;
    float* zz;
    xx = &xyz[0];
    yy = &xyz[1];
    zz = &xyz[2];

      int j = -3;
      for(int i = 0;i < NProtAtoms;i++) 
      {
      j+=3;
      // fill up the proteins q/m the host in the w position of the float4s
      protein_host[i].x = xx[j];
      protein_host[i].y = yy[j];
      protein_host[i].z = zz[j];
      }

   
    // copy the water from the host to the device
   cudaMemcpy(protein_dev, protein_host, memProtF4, cudaMemcpyHostToDevice);

}

/*

set_prot_sel_gpu_: (called from fortran as "set_prot_sel_gpu") is the fortran interface function for the 
                kernel to set an arbitrary (partial) protein structure from a selection
                to protein_sel_dev on the GPU

*/
extern "C" void set_prot_sel_gpu_(float* xyz)
{
    float* xx;
    float* yy;
    float* zz;
    xx = &xyz[0];
    yy = &xyz[1];
    zz = &xyz[2];

      int j = -3;
      for(int i = 0;i < NSelAtoms;i++) 
      {
       j+=3;
       
      // fill up the selection q/m the host in the w position of the float4s
       protein_sel_host[i].x = xx[j];
       protein_sel_host[i].y = yy[j];
       protein_sel_host[i].z = zz[j];
      }

   
    // copy the water from the host to the device
   cudaMemcpy(protein_sel_dev, protein_sel_host, memSelF4, cudaMemcpyHostToDevice);

}

extern "C" void get_shl_cos_gpu_(float *xyz, float* qrms, float* prMx, float* prMy, float* prMz, float* dCos, int* lShell,int* nShell,float* cosShell,int* binCos)
{
   int nSh;
   nSh = *nShell;

   // setup the grid to run the kernel
   //int numThreadsPerBlock = 256;
   int numThreadsPerBlock = 256;
   int numBlocks = nSh / numThreadsPerBlock + (nSh % numThreadsPerBlock == 0?0:1);
 //  printf("RUNNING ON GPU WITH %i blocks and %i threads\n",numBlocks, numThreadsPerBlock);
 //  printf("RUNNING ON GPU WITH %i Nw and %i Np \n",NwAtoms, NProtAtoms);

   float* xo;
   float* yo;
   float* zo;

   float* xh1;
   float* yh1;
   float* zh1;

   float* xh2;
   float* yh2;
   float* zh2;

   float* qo;
   float* qh1;
   float* qh2;


   int iWatShell;
   for(int i = 0;i < nSh;i++)
   {

   iWatShell = (lShell[i] * 3) - 3;

   xo = &xyz[iWatShell];
   yo = &xyz[iWatShell+1];
   zo = &xyz[iWatShell+2];
   qo = &qrms[iWatShell]; // oxygen

   xh1 = &xyz[iWatShell+3];
   yh1 = &xyz[iWatShell+4];
   zh1 = &xyz[iWatShell+5];
   qh1 = &qrms[iWatShell+3]; // h1

   xh2 = &xyz[iWatShell+6];
   yh2 = &xyz[iWatShell+7];
   zh2 = &xyz[iWatShell+8];
   qh2 = &qrms[iWatShell+6]; // h2

    water_oxy_host[i] = make_float4(*xo,*yo,*zo,*qo);
    water_h1_host[i] = make_float4(*xh1,*yh1,*zh1,*qh1);
    water_h2_host[i] = make_float4(*xh2,*yh2,*zh2,*qh2);

   }

     // copy the water from the host to the device
   cudaMemcpy(water_oxy_dev, water_oxy_host, memWaterF4, cudaMemcpyHostToDevice);
   //checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h1_dev, water_h1_host, memWaterF4, cudaMemcpyHostToDevice);
   //checkCUDAError2("cudaMemcpy(wat) -- host to device");
   cudaMemcpy(water_h2_dev, water_h2_host, memWaterF4, cudaMemcpyHostToDevice);
   //checkCUDAError2("cudaMemcpy(wat) -- host to device");

   cosShell_GPU<<<numBlocks, numThreadsPerBlock>>>(nSh,*prMx,*prMy,*prMz,*dCos,
                                                  water_oxy_dev,water_h1_dev,
                                                  water_h2_dev,
                                                  water_binfo_dev);

   // copy the p2args from the host to the device
   cudaMemcpy(water_binfo_host,water_binfo_dev,memWaterF4,cudaMemcpyDeviceToHost);
   //checkCUDAError2("cudaMemcpy -- results from device to host");

   // Sum the whole array p2arg_host array, and calculate the order parameter
    float cosWP = 0.0f;
    int ibin;
    for(int i=0;i<nSh;i++)
    {
      cosWP+=water_binfo_host[i].x;
      ibin=(int)water_binfo_host[i].y;
      binCos[i] = ibin;
    }
    *cosShell = cosWP/((float)nSh);

}

void checkCUDAError2(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}

