//////////////////////////////////////////////////////////////////
//////////////// HEADER SECTION //////////////////////////////////
#include <stdio.h>

/* This is the file pointer for handling dcd files */
FILE *dcdPtr;

void read_dcd_file_head_();
void open_dcd_(char *fileName, long int fileNameLen);
void read_dcd_(long int *nProc,long int *iAtoms, float *xyz, float *box);
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void close_dcd_(){
/*  This is a fortran wrapper for using fclose on the
    dcd file pointer. 
*/
	fclose(dcdPtr);	
}

void open_dcd_(char *fileName, long int fileNameLen){
/*  This is a fortran wrapper for using fopen on the
    binary dcd file.

    Args: char *fileName - the null terminated [eg. // CHAR(0)] character string
          long int fileNameLen - the length of the filenName string is implicity passed 
                                 as a 2nd arg.
*/
	char *filemode;
        filemode = "rb";

	/*	printf("The command to open is: %s \n",cmd); */
	dcdPtr = fopen(fileName,filemode);
	printf("   OPENED THE DCD FILE NAMED |%s| \n",fileName);	

        // now read out the full dcd file header
        read_dcd_file_head_();	
}

void read_dcd_file_head_(){
/*  This is a fortran wrapper for using pclose on the
    dcd file pointer. 
*/

  unsigned int eighty_four[1];
  unsigned int num_four[1];
  unsigned int title_size[1];
  unsigned int num_atoms[1];
  unsigned int num_frames;
  unsigned int start_frame;
  unsigned int frame_period;
  unsigned int num_fixed;
  double dTimeDelta;
  float  fTimeDelta;
  char headerbuffer[84];
  unsigned int i;
  char magicCORD[4];
  int read_val;

  // the first 4 bytes should be the number 84
  eighty_four[0] = 0;
  read_val = fread(eighty_four,sizeof(unsigned int),1,dcdPtr);
//   if (eighty_four[0] == 84)
//   {
//     printf("We have the number 84 on our hands \n");
//   }else {
//     printf("Dammit we dont even have fkin 84!!!! \n");
//   }

  // the next 4 characters should be 'CORD'
  read_val = fread(magicCORD,sizeof(char),4,dcdPtr);
   printf("   DCD MAGIC 'CORD' = %c%c%c%c \n",magicCORD[0],magicCORD[1],magicCORD[2],magicCORD[3]);

  // the next 80 characters as the header
  read_val = fread(headerbuffer,80,1,dcdPtr);

  i = *((int *) (headerbuffer + 76));
//  printf("Does this look like a CHARMM version number? = %i \n",i);

  i = *((int *) (headerbuffer + 40));
//  printf("Does this look like a we have an extra block? = %i \n",i);

  i = *((int *) (headerbuffer + 44));
//  printf("Do we have 4 dims of something? = %i \n",i);


  // the next integer should be the number of frames in the trajectory
  num_frames = *((int *) (headerbuffer));
  printf("   The number of frames in the DCD = %i \n",num_frames);

  // the next integer should be the starting frame number in the trajectory
  start_frame = *((int *) (headerbuffer + 4));
//   printf("Start frame = %i \n",start_frame);

  // the next integer should be the frame saving period in the trajectory
  frame_period = *((int *) (headerbuffer + 8));
//  printf("Frame period = %i \n",frame_period);

  // the next integer should be the number of fixed atoms
  num_fixed = *((int *) (headerbuffer + 32));
//  printf("Number of fixed atoms = %i \n",num_fixed);

  // next, get the delta value (double on xplor, float on charmm)
  fTimeDelta = *((float *)(headerbuffer + 36));
//  printf("(float) delta time = %f \n",fTimeDelta);

  // (this is just for possible xplor versions... maybe never used)
  //dTimeDelta =  *((double *)(headerbuffer + 36));
  //printf("(double) delta time = %f \n",dTimeDelta);

  // the next 4 bytes should be the number 84
  eighty_four[0] = 0;
  read_val = fread(eighty_four,sizeof(unsigned int),1,dcdPtr);
//   if (eighty_four[0] == 84)
//   {
//     printf("We have the number 84 on our hands \n");
//   }else {
//     printf("Dammit we dont even have fkin 84!!!! \n");
//   }

  // the next 4 bytes should be size of the title block
  title_size[0] = 0;
  read_val = fread(title_size,sizeof(unsigned int),1,dcdPtr);
//     printf("The total size of the title block is = %i \n", title_size[0]);

  // the next 4 bytes should be number of 80-length title blocks
  title_size[0] = 0;
  read_val = fread(title_size,sizeof(unsigned int),1,dcdPtr);
//     printf("The number of 80-length title blocks is = %i \n", title_size[0]);

  for(i=0;i<title_size[0];++i)
  {
   // the next set of 80 characters are the title 
   read_val = fread(headerbuffer,80,1,dcdPtr);
   headerbuffer[79] = '\0';
   printf("   DCD TITLE %i: %s \n",i+1,headerbuffer);
  }

  // again, the next 4 bytes should be size of the title block
  title_size[0] = 0;
  read_val = fread(title_size,sizeof(unsigned int),1,dcdPtr);
//  printf("Again, the total size of the title block should be here = %i \n", title_size[0]);

  // the next 4 bytes should be the number 84
  num_four[0] = 0;
  read_val = fread(num_four,sizeof(unsigned int),1,dcdPtr);
//   if (num_four[0] == 4)
//   {
//     printf("We have the number 4 on our hands \n");
//   }else {
//     printf("Dammit we dont even have fkin 4!!!! \n");
//   }

  // again, the next 4 bytes should be the number of atoms
  num_atoms[0] = 0;
  read_val = fread(num_atoms,sizeof(unsigned int),1,dcdPtr);
//  printf("Now, the number of atoms = %i \n", num_atoms[0]);

  // the next 4 bytes should be the number 84
  num_four[0] = 0;
  read_val = fread(num_four,sizeof(unsigned int),1,dcdPtr);
//   if (num_four[0] == 4)
//   {
//     printf("YAY, another 4! \n");
//   }else {
//     printf("Dammit we dont even have fkin 4!!!! \n");
//   }

/*

  TO DO: here, i need to write a header reading section for fixed atoms,
         but at this point i just want to have a dcd file read and 
         really dont care how i manage to scrape things together.
*/
}

void read_dcd_(long int *nProc,long int *iAtoms, float *xyz, float *box){
/*  This is a fortran wrapper for reading dcd files.  The routine will read nProc worth of 
    configurations into xyz and box lengths into box.

    Args: long int *nProc - the number of processors
          long int *iAtoms - the number of atoms in the system
          float *xyz - pointer to the "basic_atom" data type (3 reals: x,y,z)
          float *box - pointer to the "boxL" data type (3 reals: x,y,z)
*/		
	int i,j,iConf,iTot,iBox, nProc1;
	int iptr, iTot1;
	float *xp,*yp,*zp;	
	float *boxxp,*boxyp,*boxzp;
        int read_val;
        double cell_data[6];
        unsigned int forty_eight[1];
        unsigned int sizeofiAtoms[1];
        float angX,angY,angZ;
        float X_buffer[*iAtoms];
        float Y_buffer[*iAtoms];
        float Z_buffer[*iAtoms];
	boxxp = &box[0];  /* boxxp -> boxL%x */
	boxyp = &box[1];  /* boxyp -> boxL%y */ 
	boxzp = &box[2];  /* boxzp -> boxL%z */					


	iTot = (*iAtoms * 3) - 1;
	iTot1 = iTot + 1;
	iBox = 0;
	iptr = 0-iTot1;	
	nProc1 = *nProc - 1;
	for (iConf=0;iConf<= nProc1;iConf++){
	/* Set up pointers for fast pointer arithmetic using i/iBox as pointer indexes */
		iptr = iptr+iTot1;

		//xp = &xyz[iptr]; /* xp -> basic_atom%x */
		//yp = &xyz[iptr+1]; /* yp -> basic_atom%y */
		//zp = &xyz[iptr+2]; /* zp -> basic_atom%z */			

                forty_eight[0] = 0;
                read_val = fread(forty_eight,sizeof(unsigned int),1,dcdPtr);
/*                if (forty_eight[0] == 48)
                {
                  printf("We have the number 48 on our hands \n");
                }else {
                  printf("Dammit we dont even have fkin 48!!!! \n");
                }
*/
               // the next 80 characters as the header
               read_val = fread(cell_data,48,1,dcdPtr);
               
               boxxp[iBox] = (float)cell_data[0];
               boxyp[iBox] = (float)cell_data[2];
               boxzp[iBox] = (float)cell_data[5];
//               angX = (float)cell_data[1];
//               angY = (float)cell_data[3];
//               angZ = (float)cell_data[4];

//               printf("FROM DCD: boxX = %f boxY = %f bozZ = %f \n",boxxp[iBox],boxyp[iBox],boxzp[iBox]);
//               printf("FROM DCD: angX = %f angY = %f angZ = %f \n",angX,angY,angZ);
               iBox+=3;	

                forty_eight[0] = 0;
                read_val = fread(forty_eight,sizeof(unsigned int),1,dcdPtr);
/*                if (forty_eight[0] == 48)
                {
                  printf("We have the number 48 on our hands \n");
                }else {
                  printf("Dammit we dont even have fkin 48!!!! \n");
                }
*/
                 sizeofiAtoms[0] = 0;
                read_val = fread(sizeofiAtoms,sizeof(unsigned int),1,dcdPtr);
//                  printf("The size of N*sizeof(float) atoms = %i \n",sizeofiAtoms[0]);
                read_val = fread(X_buffer,sizeofiAtoms[0],1,dcdPtr);
                 sizeofiAtoms[0] = 0;
                read_val = fread(sizeofiAtoms,sizeof(unsigned int),1,dcdPtr);
//                  printf("The size of N*sizeof(float) atoms, again = %i \n",sizeofiAtoms[0]);

                 sizeofiAtoms[0] = 0;
                read_val = fread(sizeofiAtoms,sizeof(unsigned int),1,dcdPtr);
//                  printf("The size of N*sizeof(float) atoms = %i \n",sizeofiAtoms[0]);
                read_val = fread(Y_buffer,sizeofiAtoms[0],1,dcdPtr);
                 sizeofiAtoms[0] = 0;
                read_val = fread(sizeofiAtoms,sizeof(unsigned int),1,dcdPtr);
//                  printf("The size of N*sizeof(float) atoms, again = %i \n",sizeofiAtoms[0]);


                 sizeofiAtoms[0] = 0;
                read_val = fread(sizeofiAtoms,sizeof(unsigned int),1,dcdPtr);
//                  printf("The size of N*sizeof(float) atoms = %i \n",sizeofiAtoms[0]);
                read_val = fread(Z_buffer,sizeofiAtoms[0],1,dcdPtr);
                 sizeofiAtoms[0] = 0;
                read_val = fread(sizeofiAtoms,sizeof(unsigned int),1,dcdPtr);
//                  printf("The size of N*sizeof(float) atoms, again = %i \n",sizeofiAtoms[0]);

 		//read_val = fscanf(dcdPtr,"%f %f %f\n",&boxxp[iBox],&boxyp[iBox],&boxzp[iBox]);
                j=-1;
                for (i=0;i<= iTot;i+=3) {
                ++j;
             //  j = iptr + i;

//		xp = &xyz[iptr]; /* xp -> basic_atom%x */
//		yp = &xyz[iptr+1]; /* yp -> basic_atom%y */
//		zp = &xyz[iptr+2]; /* zp -> basic_atom%z */			

                 xyz[iptr + i    ] = (float)X_buffer[j];
                 xyz[iptr + i + 1] = (float)Y_buffer[j];
                 xyz[iptr + i + 2] = (float)Z_buffer[j];
//                 if (j == 0) {
//                  printf("IN DCD: x = %f y = %f z = %f \n\n", X_buffer[j], Y_buffer[j], Z_buffer[j]);
//                 }
                }
		//for (i=0;i <= iTot;i+=3){
		//	read_val = fscanf(dcdPtr,"%f %f %f",&xp[i],&yp[i],&zp[i]); 
		//}		
	}
}
