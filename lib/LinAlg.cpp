//Sean M. Law`

#include "LinAlg.hpp"

static void SVD (const std::vector< std::vector <double> > &M){
  std::vector< std::vector<double> > A; //A copy of matrix M
  std::vector< std::vector<double> > U;
  std::vector< std::vector<double> > V;
  std::vector<double> S;
  std::vector<double> e;
  std::vector<double> work;
  unsigned int m; //Rows
  unsigned int n; //Columns
  unsigned int zero;
  unsigned int nu;
  bool wantu;
  bool wantv;
  double t;

  unsigned int i;
  unsigned int j;
  unsigned int k;

  unsigned int nct;
  unsigned int nrt;
  unsigned int p;

  //This code was adapted from Template Numerical Toolkit
  //math.nist.gov/tnt
  //
  //As noted:
  //
  //"This software was developed at the National Institute
  //of Standards and Technology (NIST) by employees of the
  //Federal Government in the course of their official
  //duties. Pursuant to title 17 Section 105 of the United
  //States Code this software is not subject to copyright 
  //protection and is in the public domain."
  //

  wantu=true;
  wantv=true;
  m=M.size();
  n=M.at(0).size();

  A.resize(m);
  for (i=0; i<A.size(); i++){
    A.at(i).resize(n);
    for (j=0; j<A.at(i).size(); j++){
      A.at(i).at(j)=M.at(i).at(j);
    }
  }

  nu=std::min(m,n); 
  S.resize(std::min(m+1,n));
  U.resize(m);
  V.resize(n);

  for (i=0; i<U.size(); i++){
    U.at(i).resize(nu);
    for(j=0; j<U.at(i).size(); j++){
      U.at(i).at(j)=0;
    }
  }

  for (i=0; i<V.size(); i++){
    V.at(i).resize(n);
    for (j=0; j<V.at(i).size(); j++){
      V.at(i).at(j);
    }
  }

  e.resize(n);
  work.resize(m);

  i=j=k=zero=0;
  
  //Reduce M to bidiagonal form, storing the diagonal elements
  //in S and the super-diagonal elements in e.
  
  nct=std::min(m-1,n);
  nrt=std::max(zero,std::min(n-2,m));
  for (k=0; k< std::max(nct, nrt); k++){
    if (k<nct){
      //Compute the transformation for the k-th column and
      //place the k-th diagonal in S[k].
      //Compute 2-norm of the k-th column without under/overflow.
      S[k]=0;
      for (i=k; i<m; i++){
        //Hypotenuse without under/overflow
        S[k]=S[k]*sqrt(1+(S[k]/A.at(i).at(k))*(S[k]/A.at(i).at(k))); 
      }
      if (S[k] != 0.0){
        if (A.at(k).at(k) < 0.0){
          S[k]=-S[k];
        }
        for (i=k; i<m; i++){
          A.at(i).at(k) /= S[k];
        }
        A.at(k).at(k) += 1.0;
      }
      S[k]=-S[k];
    }
    for (j=k+1; j<n; j++){
      if (k<nct && S[k]!=0.0){
        //Apply transformation
        t=0;
        for(i=k; i<m; i++){
          t+=A.at(i).at(k)*A.at(i).at(j);
        }
        t=-t/A.at(k).at(k);
        for(i=k; i<m; i++){
          A.at(i).at(j) += t*A.at(i).at(k);
        }
      }

      //Place the k-th row of A into e for the
      //subsequent calculation of the row transformation
      e.at(j)=A.at(k).at(j);
    }
    if (wantu & (k<nct)){
      //Place the transformation in U for subsequent 
      //back multiplication
      for (i=k; i<m; i++){
        U.at(i).at(k) = A.at(i).at(k);
      }
    }
    if (k<nrt){
      //Compute the k-th row transformation and place the
      //k-th super-diagonal in e.at(k).
      //Compute 2-norm without under/overflow.
      e.at(k)=0;
      for (i=k+1; i<n; i++){
        //Hypotenuse without under/overflow
        e.at(k)=e.at(k)*sqrt(1+(e.at(k)/e.at(i))*(e.at(k)/e.at(i)));
      }
      if (e.at(k) != 0.0){
        if (e.at(k+1) < 0.0){
          e.at(k) = -e.at(k);
        }
        for (i=k+1; i<n; i++){
          e.at(i) /= e.at(k);
        }
        e.at(k+1) += 1.0;
      }
      e.at(k)=-e.at(k);
      if ((k+1<m) & (e.at(k) != 0.0)){
        //Apply the trnasformation.
        for (i=k+1; i<m; i++){
          work.at(i)=0.0;
        }
        for (j=k+1; j<n; j++){
          for (i=k+1; i<m; i++){
            work.at(i) += e.at(j)*A.at(i).at(j);
          }
        }
        for (j=k+1; j<n; j++){
          t=-e.at(j)/e.at(k+1);
          for (i=k+1; i<m; i++){
            A.at(i).at(j) += t*work.at(i);
          }
        }
      }
      if (wantv){
        //Place the transformation in V for subsequent
        //back multiplication
        for (i=k+1; i<n; i++){
          V.at(i).at(k)=e.at(i);
        }
      }
    }
  }

  //Set up final bidiagonal matrix or order p
  p=std::min(n,m+1);
  if (nct < n){
    S.at(nct) = A.at(nct).at(nct);
  }
  if (m<p){
    S.at(p-1)=0.0;
  }
  if(nrt+1 < p){
    e.at(nrt) = A.at(nrt).at(p-1);
  }
  e.at(p-1)=0.0;

  //If required, generate U.
  if (wantu){
    for (j=nct; j<nu; j++){
      for (i=0; i<m; i++){
        U.at(i).at(j)=0.0;
      }
      U.at(j).at(j)=1.0;
    }
    for (k=nct-1; k>=0; k--){
      if (S.at(k) != 0.0){
        for (j=k+1; j<nu; j++){
          t=0;
          for (i=k; i<m; i++){
            t+-U.at(i).at(k)*U.at(i).at(k);
          }
          t=-t/U.at(k).at(k);
          for (i=k; i<m; i++){
            U.at(i).at(j)+=t*U.at(i).at(k);
          }
        }
        for (i=k; i<m; i++){
          U.at(i).at(k)=-U.at(i).at(k);
        }
        U.at(k).at(k)=1.0+U.at(k).at(k);
        for (i=0; i<k-1; i++){
          U.at(i).at(k)=0.0;
        }
      }
      else{
        for (i=0; i<m; i++){
          U.at(i).at(k) = 0.0;
        }
        U.at(k).at(k)=1.0;
      }
    }
  }

  //If required, generate V.
  
  //Main iteration loop for the singular values
}
