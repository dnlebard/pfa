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
//void read_psf_file_head_(float* masses, float* charges, float* radii, bool* set_from_sel, bool* do_cluster_ana, int* res_list, int* max_atoms, int* num_cluster_res, int* atoms_sel_idx, int* totalSelected, char *selection, char *residue, long int selLenth, long int resLenth);
void read_psf_file_head_(float* masses, float* charges, float* radii, bool* set_from_sel, bool* do_cluster_ana, int* res_list, int* full_res_list, int* max_atoms, int* num_cluster_res, int* atoms_sel_idx, int* totalSelected, char *selection, char *residue, long int selLenth, long int resLenth);
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
    case 'E': radGuess = 1.73f;break;  // E0
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


//void read_psf_file_head_(float* masses, float* charges, float* radii, bool* set_from_sel, bool* do_cluster_ana, int* res_list, int* max_atoms, int* num_cluster_res, int* atoms_sel_idx, int* totalSelected, char *selection, char *residue, long int selLenth, long int resLenth){
void read_psf_file_head_(float* masses, float* charges, float* radii, bool* set_from_sel, bool* do_cluster_ana, int* res_list, int* full_res_list, int* max_atoms, int* num_cluster_res, int* atoms_sel_idx, int* totalSelected, char *selection, char *residue, long int selLenth, long int resLenth){
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
  char *resName;
  char *strCharge;
  char *strMass;
  char* delim = " ";
  int max_select;
  int nSelects;
  int nSelRes;
  int nTotalSel;
  int  nameLength;
  char** select_array = NULL;
  char** residue_array = NULL;
  int thisResID;
  int lastResID;
  int thisFullResID;
  int lastFullResID;
  bool foundAtom;
  bool foundResidue;
  int  next_res;
  int  nClusterRes;
  int  nClusterAtom;
  int  nFullRes;
  int  nFullAtom;
  int* resid;
  int* atoms_per_res;
  int* atomid;
  int* full_resid;
  int* full_atoms_per_res;
  int* full_atomid;

  //setbuf(stdout,NULL); /* set output to unbuffered */

  nFullRes      = -1;
  nClusterRes   = -1;
  next_res      = *max_atoms + 2;  // the 2 is for the resid and the atoms_per_res

  max_select = 10;
  nSelects = 0;
  nSelRes = 0;
  nTotalSel = -1;
  lastResID = -1;
  lastFullResID = -1;
//  thisResID = 0;

  if(*set_from_sel)
  {
    select_array = create_array_(selection,delim,max_select);       

    for(i = 0;select_array[i] != NULL;++i)
     {nSelects++;}

    if(*do_cluster_ana)
    {
      //printf("************************************\n");
      //printf("*** The C residue selection string: |%s| *** \n",residue);

      residue_array = create_array_(residue,delim,max_select);
      for(i = 0;residue_array[i] != NULL;++i)
      {nSelRes++;}      
      //printf("***     FOUND %i residue selections \n",nSelRes);
      //printf("************************************\n");
    }

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
// printf(" in c ::: I think the number of lines in the title is ... %i \n",num_lines);
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

//  for (i = 0;i< num_lines;++i)
//  {
//    read_val = fgets(headerbuffer,100,psfPtr);  
//    printf("REMARKSECTION: |%s| \n",headerbuffer);
//  }

  for (j = 0;j< num_lines;++j)
  { 
   read_val = fgets(headerbuffer,100,psfPtr);  
   atomLine = strtok(headerbuffer," ");

   for (i = 0;i< 7;++i)
   {
    atomLine = strtok(NULL, " ");

    if(i == 1)
    {
      thisResID = atoi(atomLine);
      thisFullResID = thisResID;
    }

    if(i == 2)
    {
     resName = atomLine;

       if(*do_cluster_ana)
       {
        foundResidue = false;
        for(k = 0;k<nSelRes;++k)
        {           
           if (strcmp(residue_array[k],resName) == 0)
           {foundResidue = true;break;}
        }

       }
    }

    if(i == 3)
     {atomName = atomLine;
      radii[j] = guess_radii_(atomName);

       if(*set_from_sel)
       {
         foundAtom = false;
         for(k = 0;k<nSelects;++k)
         {           
           if (strcmp(select_array[k],atomName) == 0)
           {
             atoms_sel_idx[++nTotalSel] = j + 1;
             foundAtom = true;break;
           }
         }

       }

      }
    if(i == 5)
     {charges[j] = atof(atomLine);}
    if(i == 6)
     {masses[j] = atof(atomLine);}
   }


/////////////////////////////
    if(*do_cluster_ana)
    {

      ///////// test for res_list ONLY ///////
      if( foundResidue && foundAtom )
      {
        if( lastResID != thisResID )
        {
          nClusterAtom = -1;
          nClusterRes++;
          lastResID = thisResID;

          //printf("nClusterAtom = %i  \n",nClusterAtom);
          //printf("nClusterRes = %i  \n",nClusterRes);
          //printf("thisResID = %i  \n",thisResID);

          resid         = &res_list[0 + (nClusterRes*next_res)];
          atoms_per_res = &res_list[1 + (nClusterRes*next_res)];
          atomid        = &res_list[2 + (nClusterRes*next_res)];

          resid[0] = thisResID;
        }

        nClusterAtom++;
        atoms_per_res[0]     = nClusterAtom + 1;
        atomid[nClusterAtom] = atoms_sel_idx[nTotalSel];
      } 
      ///////// end test for res_list ONLY ///////

      ///////// test for full_res_list ///////
      if( foundResidue )
      {
        if( lastFullResID != thisFullResID )
        {
          nFullAtom = -1;
          nFullRes++;
          lastFullResID = thisFullResID;

          //printf("nFullAtom = %i  \n",nFullAtom);
          //printf("nFullRes = %i  \n",nFullRes);
          //printf("thisFullResID = %i  \n",thisFullResID);

          full_resid         = &full_res_list[0 + (nFullRes*next_res)];
          full_atoms_per_res = &full_res_list[1 + (nFullRes*next_res)];
          full_atomid        = &full_res_list[2 + (nFullRes*next_res)];

          full_resid[0] = thisFullResID;
        }

        nFullAtom++;
        full_atoms_per_res[0]  = nFullAtom + 1;
        full_atomid[nFullAtom] = j + 1;
      } 
      ///////// end test for full_res_list ///////

    }
/////////////////////////////

  }


  if(*set_from_sel)
  {
   ++nTotalSel; // this is to account for the 0-based array

   //printf("TOTAL SELECTED IN C: %i \n",nTotalSel);
   *totalSelected = nTotalSel;

    if(*do_cluster_ana)
    {
     ++nClusterRes; // accounting for a 0-based array
     *num_cluster_res = nClusterRes;

     // make sure we free up the memory
     for(i = 0;i<nSelRes;++i)
      {free(residue_array[i]);  }      
     free(residue_array);

    }

    // make sure we free up the memory, again
    for(i = 0;i<nSelects;++i)
     {  free(select_array[i]);  }
    free(select_array);


  }

}

