#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXINPFILES 4096
#define MAXLINESIZE 4096
#define kB (1.380650424E-23)*(1.44E20) 
#define MAXGUESSBIN 10
#define TINY 1.0e-20;

//void test (double **a, int n, int *indx, double *d){
void test (){
  double *v;
  int n;
  n=15;
  v=(double *)malloc(sizeof(double));
  *v=14.0*19;
  v=NULL;
  free(v);
  v=NULL;
}

void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double vv[n+1];
  
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) {
      fprintf (stderr,"Singular matrix in routine ludcmp\n");
      exit (0);
    }
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++){
	sum -= a[i][k]*a[k][j];
      }
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
}
//#undef TINY
/* (C) Copr. 1986-92 Numerical Recipes Software *1(.|a. */

void lubksb(double **a, int n, int *indx, double *b)
{
  int i,ii=0,ip,j;
  double sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software *1(.|a. */

main (int argc, char **argv){
  int a,i,j,k,l,m,o,p,w;
  char *inp[MAXINPFILES];
  int ninp;
  FILE *file;
  char line[MAXLINESIZE];
  char fin[MAXLINESIZE],*finull=NULL;
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
  int guessbin,fast;
  time_t t1, t2;
  double **pd; //Partial derivatives
  double **am, *bm, *dm; //Matrix 
  int *indxm; //LU Decomposition 
  
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
  fast=0;
  minE=1E10;
  *fin=0;
  
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
    else if (!strcmp(argv[i],"-fast")){
      fast=1;
    }
    else if (!strcmp(argv[i],"-fin")){
      strcpy(fin,argv[++i]);
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
      fprintf(stderr,"Error: Could not allocate memory for \"data\"\n");
      exit (0);
    }
    if ((ebw=(double ***)malloc((ninp)*sizeof(double **)))==NULL){
      fprintf(stderr,"Error: Could not allocate memory for \"ebw\"\n");
      exit (0);
    }
    if ((n=(int *)malloc((ninp)*sizeof (int)))==NULL){
      fprintf(stderr,"Error: Could not allocate memory for \"n\"\n");
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
	  fprintf(stderr,"Error: Could not allocate memory for \"xi\"\n");
	  exit (0);
	}
	if ((ebw[i]=(double **)malloc(l*sizeof(double *)))==NULL){
	  fprintf (stderr,"Error: Could not allocate memory for \"ebw\"\n");
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
	    fprintf (stderr,"Error: Could not allocate memory for \"ebw\"\n");
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
    fprintf (stderr,"Error: Could not allocate memory for \"ebf\"\n");
    exit (0);
  }
  if ((ebf2=(double *)malloc((ninp)*sizeof(double)))==NULL){
    fprintf (stderr,"Error: Could not allocate memory for \"ebf2\"\n");
    exit (0);
  }
  if ((f=(double *)malloc((ninp)*sizeof(double)))==NULL){
    fprintf (stderr,"Error: Could not allocate memory for \"f\"\n");
    exit (0);
  }
  if ((fold=(double *)malloc((ninp)*sizeof(double)))==NULL){
    fprintf (stderr,"Error: Could not allocate memory for \"fold\"\n");
    exit (0);
  }
  if ((ftmp=(double *)malloc((ninp)*sizeof(double)))==NULL){
    fprintf (stderr,"Error: Could not allocate memory for \"ftmp\"\n");
    exit (0);
  }
  

  for (i=0; i<ninp; i++){
    f[i]=0.0;
    fold[i]=0.0;
    ftmp[i]=0.0;
  }
  if (*fin != 0){
    //Use input f values
    fprintf (stderr, "Using -fin %s\n",fin);
    if (verbose) printf ("#Using -fin %s\n",fin);
    finull=strtok(fin,":=");
    for(i=0; i<ninp && finull != NULL; i++){
      fprintf (stderr, "Setting f[%d] = %f\n",i,atof(finull));
      if (verbose) printf ("#Setting f[%d] = %f\n",i,atof(finull));
      f[i]=atof(finull);
      fold[i]=atof(finull);
      ftmp[i]=atof(finull);
      finull=strtok(NULL,":=");
    }
  }
  
  //*****************************
  //Pre-guess f[i] with less windows
  /*
  if (all || guess){
    o=0;
    while (o < ninp && o != ninp-1){
      p=o+guessbin-1;
      if (o >= ninp-guessbin){
	while (o > ninp-guessbin){
	  o--;
	}
	p=ninp-1;
      }
      
    repeat:

      for (m=1; m<=iter; m++){
	if (m != 0 && m % 10000 == 0){
	  fprintf (stderr,"#Pre-Guess Iteration %d: f[%d] --> f[%d]\n",m,o,p);
	  if (verbose) printf ("#Pre-guess Iteration %d: f[%d] --> f[%d]\n",m,o,p);
	}
	for (k=o; k<=p; k++){
	  ebfk=0.0;
	  for (i=o; i<=p; i++){
	    for(l=0; l<n[i]; l++){
	      bottom=0.0;
	      for (j=o; j<=p; j++){
		bottom=bottom+ebw[i][l][j]*(exp(ftmp[j]/(kB*T)))*n[j];
	      }
	      ebfk=ebfk+ebw[i][l][k]/bottom;
	    }
	  }
	  
	  ftmp[k]=-kB*T*log(ebfk);
	}
	
	for(k=o+1; k<=p; k++){
	  ftmp[k]-=ftmp[o];
	}
	
	conv=1;
	for (k=o; k<=p; k++){
	  delta=fabs(ftmp[k]-fold[k]);
	  if (delta>=tol) conv=0;
	  fold[k]=ftmp[k];
	}
	if (conv>0) break;
	fflush(stdout);
      }
      for (k=o; k<=p; k++){
	//Check and see if ftmp==NaN
	//If yes, expand the guessbin
	if (ftmp[k] != ftmp[k]){
	  for (i=o; i<=p; i++){
	    if (ftmp[k] != ftmp[k]){
	      ftmp[k]=0.0; //Reset ftmp[k] in guessbin range
	    }
	  }
	  if ( o>0 ) o--;
	  fprintf (stderr,"#Expanding guessbin: f[%d] --> f[%d]\n",o,p);
	  if (verbose) printf ("#Expanding guessbin: f[%d] --> f[%d]\n",o,p);
	  if (p-o == MAXGUESSBIN){
	    fprintf (stderr,"Error: Guessbin size has exceeded the limit\n");
	    if (verbose) printf ("Error: Guessbin size has exceeded the limit\n");
	  }
	  goto repeat;
	}
      }
    
      for (k=o; k<=p; k++){
	if (o == 0){
	  f[k]=ftmp[k];
	}
	else{
	  f[k]=(f[o]-ftmp[o])+ftmp[k];
	}
	fprintf (stderr, "#Pre-Guess: f[%d] = %f\n",k, f[k]);
	if (verbose) printf ("#Pre-Guess: f[%d] = %f\n",k, f[k]);
      }
      
      o=p;
    }
  }
  */
  //***********


  if (fast){
    if((am=(double **)malloc((ninp+1)*sizeof(double *)))==NULL){//Correct!
      fprintf (stderr,"Error: Could not allocate memory for \"am\"\n");
      exit (0);
    }
    for (i=0; i<=ninp; i++){
      if((am[i]=(double *)malloc((ninp+1)*sizeof(double)))==NULL){
	fprintf (stderr,"Error: Could not allocate memory for \"am\"\n");
	exit (0);
      }
    }
    if((bm=(double *)malloc((ninp+1)*sizeof(double)))==NULL){
      fprintf (stderr,"Error: Could not allocate memory for \"bm\"\n");
      exit (0);
    }
    if((indxm=(int *)malloc((ninp+1)*sizeof(int)))==NULL){
      fprintf (stderr,"Error: Could not allocate memory for \"indxm\"\n");
      exit (0);
    }
    if((dm=(double *)malloc(sizeof(double)))==NULL){
      fprintf (stderr,"Error: Could not allocate memory for \"dm\"\n");
      exit (0);
    }
  }

  //*****************
  //Add one more bin at a time
  if (all || guess){
    o=0;
    p=1;
    while (p <= ninp-1){
      for (m=1; m<=iter; m++){
	if (m != 0 && m % 10000 == 0){
	  fprintf (stderr,"#Pre-Guess Iteration %d: f[%d] --> f[%d]\n",m,o,p);
	  if (verbose) printf ("#Pre-guess Iteration %d: f[%d] --> f[%d]\n",m,o,p);
	}
	for (k=o; k<=p; k++){
	  ebfk=0.0;
	  for (i=o; i<=p; i++){
	    for(l=0; l<n[i]; l++){
	      bottom=0.0;
	      for (j=o; j<=p; j++){
		bottom=bottom+ebw[i][l][j]*(exp(ftmp[j]/(kB*T)))*n[j];
	      }
	      ebfk=ebfk+ebw[i][l][k]/bottom;
	    }
	  }
	  
	  ftmp[k]=-kB*T*log(ebfk);
	}
	
	for(k=o+1; k<=p; k++){
	  ftmp[k]-=ftmp[o];
	}
	
	conv=1;
	for (k=o; k<=p; k++){
	  delta=fabs(ftmp[k]-fold[k]);
	  if (delta>=tol) conv=0;
	  fold[k]=ftmp[k];
	}
	if (conv>0) break;
	fflush(stdout);
	//Accelerate Convergence
	
      }
      for (k=o; k<=p; k++){
	//Check and see if ftmp==NaN
	if (ftmp[k] != ftmp[k]){
	  for (i=o; i<=p; i++){
	    if (ftmp[k] != ftmp[k]){
	      ftmp[k]=0.0; //Reset ftmp[k] in guessbin range
	    }
	  }
	}
      }
    
      for (k=o; k<=p; k++){
	if (o == 0){
	  f[k]=ftmp[k];
	}
	else{
	  f[k]=(f[o]-ftmp[o])+ftmp[k];
	}
	fprintf (stderr, "#Pre-Guess: f[%d] = %f\n",k, f[k]);
	if (verbose) printf ("#Pre-Guess: f[%d] = %f\n",k, f[k]);
      }
      fprintf (stderr, "#%d Iterations\n\n", m);
      if (verbose) printf ("#%d Iterations\n\n",m);

      p++;
    }
  }

  //***********

  if (all){
    for (i=0; i<ninp; i++){
      if (f[i] != f[i]) f[i]=0.0; //Means that f[i]=NaN
    }
  }
  
  if (fast){
    //Set up partial derivatives
    if ((pd=(double **)malloc((ninp)*sizeof(double *)))==NULL){
      fprintf(stderr,"Error: Could not allocate memory for \"pd\"\n");
      exit (0);
    }
    for (i=0; i<ninp; i++){
      if ((pd[i]=(double *)malloc(sizeof(double)))==NULL){
	      fprintf(stderr,"Error: Could not allocate memory for \"pd\"\n");
	      exit (0);
      }
    }
  }
  

  if (all || noguess){
    for (m=1; m<=iter; m++){
      fprintf (stderr,"Iteration %d,", m);
      if (verbose) printf ("#Iteration %d,", m);
      t1=time(NULL);
      for (k=0; k<ninp; k++){
	ebfk=0.0;
	for (i=0; i<ninp; i++){
	  for(l=0; l<n[i]; l++){
	    bottom=0.0;
	    for (j=0; j<ninp; j++){
	      //bottom=bottom+ebw[i][l][j]*(exp(f[j]/(kB*T)))*n[j];
	      bottom=bottom+ebw[i][l][j]*(exp(fold[j]/(kB*T)))*n[j];
	    }
	    ebfk=ebfk+ebw[i][l][k]/bottom;
	  }
	}
	
	f[k]=-kB*T*log(ebfk);
      }
      
      for(k=1; k<ninp; k++){
	f[k]-=f[0];
      }
      //Accelerate convergence
      if (fast){
	//Get Partial Derivatives
	for (k=0; k<ninp; k++){ //Take derivative of Fk WRT
	  for (a=0; a<ninp; a++){ //Fa
	    ebfk=0.0;
	    for (w=0; w<ninp; w++){
	      for(l=0; l<n[w]; l++){
		bottom=0.0;
		for (j=0; j<ninp; j++){
		  bottom=bottom+
		    (ebw[w][l][j]*(exp(fold[j]/(kB*T)))*n[j]);
		}
		ebfk=ebfk+(ebw[w][l][a]*ebw[w][l][k])/(bottom*bottom);
	      }
	    }
	    pd[k][a]=n[a]*(exp((f[k]+fold[a])/(kB*T)))*ebfk;
	  }
	}
	
	//Generate Matrix
	for (k=0; k<ninp; k++){ //For each Fk, 
	  //for (i=0; i<ninp; i++){
	  bm[k+1]=f[k];
	  for (a=0; a<ninp; a++){
	    bm[k+1]-=pd[k][a]*fold[a];
	    if (k != a){
	      am[k+1][a+1]=-pd[k][a];
	    }
	    else{
	      am[k+1][a+1]=1-pd[k][a];
	    }
	  }
	    //}
	}
	
	//Calculate Solution for the desired F values
	ludcmp(am, ninp, indxm, dm);
	lubksb(am, ninp,indxm, bm);
	for(k=1; k<ninp; k++){
	  bm[k+1]-=bm[1];
	}
	for (k=0; k<ninp; k++){
	  fprintf (stderr, "NEW[%d] = %f, OLD[%d] = %f\n", k+1, bm[k+1], k, f[k]);
	  f[k]=bm[k+1];
	}
      }
      
      conv=1;
      for (k=0; k<ninp; k++){
	      delta=fabs(f[k]-fold[k]);
	      if (delta>=tol) conv=0;
	      fold[k]=f[k];
      }
      t2=time(NULL);
      if (difftime(t2,t1) > 60){
	      fprintf (stderr, " %f minutes\n", difftime(t2,t1)/60);
	      if (verbose) printf (" %f minutes\n", difftime(t2,t1)/60);
      }
      else{
	      fprintf (stderr, " %f seconds\n", difftime(t2,t1));
	      if (verbose) printf (" %f seconds\n", difftime(t2,t1));
      }
      
      if (conv>0) break;
      fflush(stdout);
     
    }
    for (k=0; k<ninp; k++){
      fprintf (stderr, "f[%d] = %f\n", k, f[k]);
      if (verbose) printf ("#f[%d] = %f\n", k, f[k]);
    }

    //Out put f values for fin option
    fprintf (stderr, "-fin %f", f[0]);
    if (verbose) printf ("#-fin %f",f[0]);
    for (k=1; k<ninp; k++){
      fprintf (stderr, ":%f", f[k]);
      if (verbose) printf (":%f",f[k]);
    }
    fprintf (stderr, " \n");
    if (verbose) printf (" \n");
  }
  
  //Get binsize, bin xi for each window
  binsize=(max-min)/bins;
  if ((prob=(double *)malloc((bins)*sizeof(double)))==NULL){
    fprintf (stderr,"Error: Could not allocate memory for \"prob\"\n");
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
    fprintf (stderr,"Error: Could not allocate memory for \"freeE\"\n");
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
