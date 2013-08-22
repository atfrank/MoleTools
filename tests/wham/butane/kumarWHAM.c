#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096
#define kB (1.380650424E-23)*(1.44E20) 
#define MAXGUESSBIN 10

main (int argc, char **argv){
  int i,j,k,l,m,o,p;
  char *inp[MAXINPFILES];
  int ninp;
  FILE *file;
  char line[MAXLINESIZE];
  char *word;
  double min,max;
  double **xi,***ebw;
  double T;
  int bins,b,*n;
  double binsize;
  double *prob, probability;
  double normal, *freeE, minE;
  int iter;
  double *ebf,*ebf2,*f,*fold, *ftmp;
  double ebfk,bottom,delta,tol;
  int conv;
  int verbose,guess,all,noguess;
  int guessbin;
  time_t t1, t2;
  
  ninp=0;
  min=1E10;
  max=-1E10;
  T=300.0;
  iter=1000000;
  tol=0.00001;
  bins=10;
  normal=0.0;
  verbose=0;
  guess=0;
  noguess=0;
  all=1;
  guessbin=4;
  minE=1E10;
  
  for (i=1; i<argc; i++){
    if (!strcmp(argv[i],"-help")){

    }
    else if (!strcmp(argv[i],"-option")){

    }
    else if (!strcmp(argv[i],"-iter")){
      iter=atoi(argv[++i]); 
    }
    else if (!strcmp(argv[i],"-tol")){
      tol=atof(argv[++i]);
    }
    else if (!strcmp(argv[i],"-bins")){
      bins=atoi(argv[++i]);
    }
    else if (!strcmp(argv[i],"-guess")){
      guess=1;
      all=0;
      noguess=0;
    }
    else if (!strcmp(argv[i],"-guessbin")){
      guessbin=atoi(argv[++i]);
    }
    else if (!strcmp(argv[i],"-all")){
      all=1;
      guess=1;
      noguess=1;
    }
    else if (!strcmp(argv[i],"-noguess")){
      all=0;
      guess=0;
      noguess=1;
    }
    else if (!strcmp(argv[i],"-verbose")){
      verbose=1;
    }
    else{
      inp[ninp]=(char *)malloc(strlen(argv[i])+1);
      strcpy(inp[ninp++],argv[i]);
    }
  }
  
  if (ninp<=0){
    fprintf(stderr,"Please provide input file\n");
    exit(1);
  }
  else { //Read input files
    if ((xi=(double **)malloc((ninp)*sizeof(double *)))==NULL){
      printf("Error: Could not allocate memory for \"data\"\n");
      exit (0);
    }
    if ((ebw=(double ***)malloc((ninp)*sizeof(double **)))==NULL){
      printf("Error: Could not allocate memory for \"ebw\"\n");
      exit (0);
    }
    if ((n=(int *)malloc((ninp)*sizeof (int)))==NULL){
      printf("Error: Could not allocate memory for \"n\"\n");
      exit (0);
    }

    for (i=0; i<ninp; i++){
      file=fopen(inp[i],"r");
      if (file !=NULL){
        l=0;
        while (fgets(line,MAXLINESIZE,file)!=NULL){
	  if (memchr(line,'#',1)== NULL && memchr(line,'\n',1)==NULL){
	    l++;
	  }
        }
        n[i]=l;
	if ((xi[i]=(double *)malloc(l*sizeof(double)))==NULL){
	  printf("Error: Could not allocate memory for \"xi\"\n");
	  exit (0);
	}
	if ((ebw[i]=(double **)malloc(l*sizeof(double *)))==NULL){
	  printf ("Error: Could not allocate memory for \"ebw\"\n");
	  exit (0);
	}
        fclose(file);
      }
      else{
        perror(inp[i]); //Error opening file
	exit(0);
      }
    }
    
    for (i=0; i<ninp; i++){
      file=fopen(inp[i],"r");
      l=-1;
      while (fgets(line,MAXLINESIZE, file) !=NULL){
	if (memchr(line,'#',1)== NULL && memchr(line,'\n',1)==NULL){
	  l++;
	  if ((ebw[i][l]=(double *)malloc((ninp)*sizeof(double)))==NULL){
	    printf ("Error: Could not allocate memory for \"ebw\"\n");
	    exit (0);
	  }
	  
	  j=0; //Count split
	  word=strtok(line,"\t \n");
	  while (word != NULL){
	    j++;
	    if (j > 2){//This is correct!
	      if (j-3 == ninp){//This is correct!
		//Number of windows must not exceed number of input files
		break;
	      }
	      ebw[i][l][j-3]=exp(-(atof(word))/(kB*T));
	      //printf ("%f %d %d %d\n", ebw[i][l][j-3], i,l,j-3);
	    }
	    else if (j==2){
	      xi[i][l]=atof(word);
	      if (xi[i][l]<min) min=xi[i][l];
	      if (xi[i][l]>max) max=xi[i][l];
	    }
	    word=strtok(NULL,"\t \n");
	  }
	}
      }      
      fclose(file);
    }
  }
  
  if ((ebf=(double *)malloc((ninp)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"ebf\"\n");
    exit (0);
  }
  if ((ebf2=(double *)malloc((ninp)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"ebf2\"\n");
    exit (0);
  }
  if ((f=(double *)malloc((ninp)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"fact\"\n");
    exit (0);
  }
  if ((fold=(double *)malloc((ninp)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"fact\"\n");
    exit (0);
  }
  if ((ftmp=(double *)malloc((ninp)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"ftmp\"\n");
    exit (0);
  }
  
  for (i=0; i<ninp; i++){
    f[i]=0.0;
    fold[i]=0.0;
    ftmp[i]=0.0;
  }

  if (all || noguess){
    for (m=1; m<=iter; m++){
      fprintf (stderr,"Iteration %d,", m);
      if (verbose) printf ("#Iteration %d,", m);
      for (k=0; k<ninp; k++){
	for (i=0; ;<ninp; i++){
	  
	}
	exp(-fj)
      }
    }
      
      
    for (k=0; k<ninp; k++){
      fprintf (stderr, "f[%d] = %f\n", k, f[k]);
      if (verbose) printf ("#f[%d] = %f\n", k, f[k]);;
    }
  }
  
  //Get binsize, bin xi for each window
  binsize=(max-min)/bins;
  if ((prob=(double *)malloc((bins)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"prob\"\n");
    exit (0);
  }
  for(b=0; b<bins; b++){
    prob[b]=0.0;
  }

  for (i=0; i<ninp; i++){
    for (l=0; l<n[i]; l++){
      b=0;
      while (xi[i][l] > min+binsize*(b+1)){
	b++;
      }
      probability=0.0;
      for (j=0; j<ninp; j++){
	probability+=n[j]*ebw[i][l][j]*exp(f[j]/(kB*T));
      }
      probability=(1.0/probability);
      //printf ("%f\n", probability);
      prob[b]+=probability;
    }
  }

  for (b=0; b<bins; b++){
    normal+=prob[b];
  }

  if ((freeE=(double *)malloc((bins)*sizeof(double)))==NULL){
    printf ("Error: Could not allocate memory for \"freeE\"\n");
    exit (0);
  }
  for (b=0; b<bins; b++){
    prob[b]=prob[b]/normal;
    freeE[b]=-kB*T*log(prob[b]);
    if (freeE[b] < minE) minE=freeE[b];
  }
  
  printf ("\n");

  for (b=0; b<bins; b++){
    freeE[b]=freeE[b]-minE;
    printf ("%f %f %f\n",min+b*binsize+binsize/2.0, freeE[b], prob[b]);
  }

}
