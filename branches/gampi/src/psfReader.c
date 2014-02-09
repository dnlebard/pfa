//////////////////////////////////////////////////////////////////
//////////////// HEADER SECTION //////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* This is the file pointer for handling psf files */
FILE *psfPtr;

void close_psf_();
void open_psf_(char *fileName, long int fileNameLen);
void read_psf_file_head_(float* masses, float* charges, float* radii, bool* set_from_sel, int* atoms_sel_idx, int* totalSelected, char *selection, long int selLenth);
char** create_array_(char* anystring, char* delimiters, int max_length);
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void close_psf_(){
/*  This is a fortran wrapper for using fclose on the
    psf file pointer. 
*/
	fclose(psfPtr);	
}

void open_psf_(char *fileName, long int fileNameLen){
/*  This is a fortran wrapper for using fopen on the
    binary psf file.

    Args: char *fileName - the null terminated [eg. // CHAR(0)] character string
          long int fileNameLen - the length of the filenName string is implicity passed 
                                 as a 2nd arg.
*/
	char *filemode;
        filemode = "r";

	/*	printf("The command to open is: %s \n",cmd); */
	psfPtr = fopen(fileName,filemode);
        //printf("Opened the psf file named |%s| \n",fileName);	

}

float guess_radii_(char *atomName) {
/* this is a helper function passes back a good guess at the atomic radius
   based on the atom name from the psf file.  hopefully, a real periodic 
   table or topology file read will replace this call.
*/

  char *normalElements = "HCNSOP";
  char atomLetter;
  float radGuess; 
  atomLetter = atomName[0];

//  printf("character 1 in atomName is ... |%c| \n",atomLetter);


  switch (atomLetter) {
    case 'C': radGuess = 1.7f;break;
    case 'N': radGuess = 1.55f;break;
    case 'S': radGuess = 1.8f;break;
    case 'O': radGuess = 1.52f;break;
    case 'P': radGuess = 1.80f;break;
    case 'F': radGuess = 1.47f;break;
    case 'H': radGuess = 1.2f;break;
    case 'M': radGuess = 1.73f;break;
    default: printf("CANT FIND ATOM NAMED |%c|\n",atomLetter);radGuess = 1.7;
  }

//  printf("radGuess = |%f| \n",radGuess);

  return radGuess;

}

char** create_array_(char* anystring, char* delimiters, int max_length) {
/*
  this is a helper function that will return an array of strings (aka, a double character pointer)
  based on a given delimiter
*/
  char** tokens = NULL;
  char* workspace = NULL;
  char* token = NULL;

  int i = 0;

//  printf("in create_array()... anystring = |%s| \n",anystring);
  tokens = malloc(sizeof(char*)*max_length);
  if(tokens == NULL)
   return NULL;

  workspace = malloc(sizeof(char) * strlen(anystring) + 1);
  if(workspace == NULL)
   return NULL;

  strcpy(workspace,anystring);
  for(i=0;i<max_length;++i)
   { tokens[i] = NULL; }

  token = strtok(workspace,delimiters);
  i = 0;

//  printf("ABOUT TO LOOP OVER TOKENS IN STRING \n");
//  printf("  INITIAL TOKEN = |%s| \n",token);

 while((i < (max_length - 1)) && (token != NULL))
 {
   tokens[i] = malloc(sizeof(char) * strlen(token) + 1);
   if(tokens[i] != NULL)
   {
     strcpy(tokens[i],token);
//     printf("     tokens[%i] = |%s| \n",i,tokens[i]);
     i++;
     token = strtok(NULL,delimiters);
   }

 }

free(workspace);
return tokens;
}


void read_psf_file_head_(float* masses, float* charges, float* radii, bool* set_from_sel, int* atoms_sel_idx, int* totalSelected, char *selection, long int selLenth){
/*  This is a fortran wrapper for reading the header of the
    psf file pointer. 
*/


  char headerbuffer[100];
  unsigned int i, j, k;
  char *read_val;
  int num_lines;
  char *token;
  char *titleLine;
  char *atomLine;
  char *atomName;
  char *strCharge;
  char *strMass;
  char* delim = " ";
  int max_select;
  int nSelects;
  int nTotalSel;
  int  nameLength;
  char** select_array = NULL;


  max_select = 10;
  nSelects = 0;
  nTotalSel = -1;
  if(*set_from_sel)
  {
//  printf("************************************\n");
//  printf("*** The C selection string: |%s| *** \n",selection);
    select_array = create_array_(selection,delim,max_select);   
//  printf("***     FOUND %i selections \n",nSelects);
    for(i = 0;select_array[i] != NULL;++i)
     {nSelects++;}
//  printf("************************************\n");
   }

  //float charge;
  //float mass;
  //float radius;

  // the first 4 bytes should be the number 84
  //headerbuffer[100] = '\0';
  read_val = fgets(headerbuffer,100,psfPtr);  
  //printf("LINE 1: %s \n",headerbuffer);
  read_val = fgets(headerbuffer,100,psfPtr);  
  //printf("LINE 2: %s \n",headerbuffer);
  read_val = fgets(headerbuffer,100,psfPtr);  
  //printf("LINE 3: |%s| \n",headerbuffer);

 titleLine = strtok(headerbuffer," ");
 num_lines = atoi( titleLine );
// printf("I think the number of lines in the title is ... %i \n",num_lines);
// while(titleLine !=NULL)
// {
//  printf ("%s\n",titleLine);
//  titleLine = strtok(NULL, " ");
// }

  for (i = 0;i< num_lines;++i)
  {
    read_val = fgets(headerbuffer,100,psfPtr);  
    //printf("REMARKSECTION: |%s| \n",headerbuffer);
  }

    read_val = fgets(headerbuffer,100,psfPtr);  
    //printf("POSTREMARK: |%s| \n",headerbuffer);

    read_val = fgets(headerbuffer,100,psfPtr);  
    //printf("POSTREMARK: |%s| \n",headerbuffer);

   titleLine = strtok(headerbuffer," ");
   num_lines = atoi( titleLine );
   //printf("I think the number of lines in the atomsection is ... %i \n",num_lines);

//  for (i = 0;i< num_lines;++i)
//  {
//    read_val = fgets(headerbuffer,100,psfPtr);  
    //printf("REMARKSECTION: |%s| \n",headerbuffer);
//  }

  for (j = 0;j< num_lines;++j)
  { 
   read_val = fgets(headerbuffer,100,psfPtr);  
   atomLine = strtok(headerbuffer," ");
   for (i = 0;i< 7;++i)
   {
    atomLine = strtok(NULL, " ");
    if(i == 3)
     {atomName = atomLine;
      radii[j] = guess_radii_(atomName);

       if(*set_from_sel)
       {
         for(k = 0;k<nSelects;++k)
         {
           if (strcmp(select_array[k],atomName) == 0)
           {atoms_sel_idx[++nTotalSel] = j + 1;break;}
         }

       }

      }
    if(i == 5)
     {charges[j] = atof(atomLine);}
    if(i == 6)
     {masses[j] = atof(atomLine);}
   }
  }


  if(*set_from_sel)
  {
   ++nTotalSel; // this is to account for the 0-based array
   //printf("TOTAL SELECTED IN C: %i \n",nTotalSel);
   *totalSelected = nTotalSel;

  // make sure we free up the memory
  for(i = 0;i<nSelects;++i)
   {  free(select_array[i]);  }
  free(select_array);
  }

}

