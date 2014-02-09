#include <stdio.h>

/* This is the file pointer for handling gzipped files */
FILE *gzPtr;


void close_gz_(){
/*  This is a fortran wrapper for using pclose on the
    gzipped file pointer. 
*/
	pclose(gzPtr);	
}

void open_gz_(char *fileName, long int fileNameLen){
/*  This is a fortran wrapper for using popen on the
    gzipped file pointer.

    Args: char *fileName - the null terminated [eg. // CHAR(0)] character string
          long int fileNameLen - the length of the filenName string is implicity passed 
                                 as a 2nd arg.
*/
	char cmd[100];

	sprintf(cmd,"zcat %s",fileName);

	/*	printf("The command to open is: %s \n",cmd); */
	gzPtr = popen(cmd,"r");
	printf("Opened the GZMDCRD file named |%s| \n",fileName);		
}

void read_gz_(long int *nProc,long int *iAtoms, float *xyz, float *box){
/*  This is a fortran wrapper for using fscanf on the
    gzipped file pointer.  The routine will read nProc worth of 
    configurations into xyz and box lengths into box.

    Args: long int *nProc - the number of processors
          long int *iAtoms - the number of atoms in the system
          float *xyz - pointer to the "basic_atom" data type (3 reals: x,y,z)
          float *box - pointer to the "boxL" data type (3 reals: x,y,z)
*/		
	int i,iConf,iTot,iBox, nProc1;
	int iptr, iTot1;
        int read_val;
	float *xp,*yp,*zp;	
	float *boxxp,*boxyp,*boxzp;

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

		xp = &xyz[iptr]; /* xp -> basic_atom%x */
		yp = &xyz[iptr+1]; /* yp -> basic_atom%y */
		zp = &xyz[iptr+2]; /* zp -> basic_atom%z */			

		for (i=0;i <= iTot;i+=3){
			read_val = fscanf(gzPtr,"%f %f %f",&xp[i],&yp[i],&zp[i]); 
		}		
		read_val = fscanf(gzPtr,"%f %f %f\n",&boxxp[iBox],&boxyp[iBox],&boxzp[iBox]);
		iBox+=3;	
	}
}
