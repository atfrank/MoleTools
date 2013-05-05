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

  unsigned int pp;
  unsigned iter;
  unsigned int kase;
  double eps;

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
      //place the k-th diagonal in S.at(k).
      //Compute 2-norm of the k-th column without under/overflow.
      S.at(k)=0;
      for (i=k; i<m; i++){
        //S[k]=hypot(S[k],A[i][k]);
        //Hypotenuse without under/overflow
        S.at(k)=S.at(k)*sqrt(1+(A.at(i).at(k)/S.at(k))*(A.at(i).at(k)/S.at(k))); 
      }
      if (S.at(k) != 0.0){
        if (A.at(k).at(k) < 0.0){
          S.at(k)=-S.at(k);
        }
        for (i=k; i<m; i++){
          A.at(i).at(k) /= S.at(k);
        }
        A.at(k).at(k) += 1.0;
      }
      S.at(k)=-S.at(k);
    }
    for (j=k+1; j<n; j++){
      if (k<nct && S.at(k)!=0.0){
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
        e.at(k)=e.at(k)*sqrt(1+(e.at(i)/e.at(k))*(e.at(i)/e.at(k)));
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
            t+=U.at(i).at(k)*U.at(i).at(k);
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
  
  if (wantv){
    for (k=n-1; k>=0; k--){
      if ((k>nrt) & (e.at(k) != 0.0)){
        for (j=k+1; j<nu; j++){
          t=0;
          for(i=k+1; i<n; i++){
            t+=V.at(i).at(k)*V.at(i).at(j);
          }
          t=-t/V.at(k+1).at(k);
          for(i=k+1; i<n; i++){
            V.at(i).at(j)+=t*V.at(i).at(k);
          }
        }
      }
      for (i=0; i<n; i++){
        V.at(i).at(k)=0.0;
      }
      V.at(k).at(k)=1.0;
    }
  }
  
  //Main iteration loop for the singular values

  pp=p-1;
  iter=0;
  eps=pow(2.0,-52.0);

  unsigned int negativeOne;
  negativeOne=-1;

  while (p > 0) {
    k=0;
    kase=0;

    // Here is where a test for too many iterations would go.

    // This section of the program inspects for
    // negligible elements in the s and e arrays.  On
    // completion the variables kase and k are set as follows.

    // kase = 1     if s(p) and e.at(k-1) are negligible and k<p
    // kase = 2     if s(k) is negligible and k<p
    // kase = 3     if e.at(k-1) is negligible, k<p, and
    //              s(k), ..., s(p) are not negligible (qr step).
    // kase = 4     if e(p-1) is negligible (convergence).

    for (k = p-2; k >= negativeOne; k--) {
      if (k == negativeOne) {
        break;
      }
      if (abs(e.at(k)) <= eps*(abs(S.at(k)) + abs(S.at(k+1)))) {
        e.at(k) = 0.0;
        break;
      }
    }
    if (k == p-2) {
      kase = 4;
    } 
    else {
      unsigned int ks;
      for (ks = p-1; ks >= k; ks--) {
        if (ks == k) {
          break;
        }
        t = (ks != p ? abs(e.at(ks)) : 0.) + (ks != k+1 ? abs(e.at(ks-1)) : 0.);
        if (abs(S.at(ks)) <= eps*t)  {
          S.at(ks) = 0.0;
          break;
        }
      }
      if (ks == k) {
        kase = 3;
      } 
      else if (ks == p-1) {
        kase = 1;
      }
      else {
        kase = 2;
        k = ks;
      }
    }
    k++;
   
    // Perform the task indicated by kase.
 
    switch (kase) {
 
      // Deflate negligible s(p).
 
      case 1: {
        double f = e.at(p-2);
        e.at(p-2) = 0.0;
        for (j = p-2; j >= k; j--) {
          //t = hypot(S.at(j),f);
          //Hypotenuse without under/overflow
          t = S.at(j)*sqrt(1+(f/S.at(j))*(f/S.at(j)));
          double cs = S.at(j)/t;
          double sn = f/t;
          S.at(j) = t;
          if (j != k) {
            f = -sn*e.at(j-1);
            e.at(j-1) = cs*e.at(j-1);
          }
          if (wantv) {
            for (i = 0; i < n; i++) {
              t = cs*V.at(i).at(j) + sn*V.at(i).at(p-1);
              V.at(i).at(p-1) = -sn*V.at(i).at(j) + cs*V.at(i).at(p-1);
              V.at(i).at(j) = t;
            }
          }
        }
      }
      break;
 
      // Split at negligible s(k).
 
      case 2: {
        double f = e.at(k-1);
        e.at(k-1) = 0.0;
        for (j = k; j < p; j++) {
          //t = hypot(S.at(j),f);
          //Hypotenuse without under/overflow
          t = S.at(j)*sqrt(1+(f/S.at(j))*(f/S.at(j)));
          double cs = S.at(j)/t;
          double sn = f/t;
          S.at(j) = t;
          f = -sn*e.at(j);
          e.at(j) = cs*e.at(j);
          if (wantu) {
            for (i = 0; i < m; i++) {
              t = cs*U.at(i).at(j) + sn*U.at(i).at(k-1);
              U.at(i).at(k-1) = -sn*U.at(i).at(j) + cs*U.at(i).at(k-1);
              U.at(i).at(j) = t;
            }
          }
        }
      }
      break;
 
      // Perform one qr step.
 
      case 3: {
        // Calculate the shift.
        double scale = std::max(std::max(std::max(std::max(
          abs(S.at(p-1)),abs(S.at(p-2))),abs(e.at(p-2))), 
          abs(S.at(k))),abs(e.at(k)));
        double sp = S.at(p-1)/scale;
        double spm1 = S.at(p-2)/scale;
        double epm1 = e.at(p-2)/scale;
        double sk = S.at(k)/scale;
        double ek = e.at(k)/scale;
        double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
        double c = (sp*epm1)*(sp*epm1);
        double shift = 0.0;
        if ((b != 0.0) | (c != 0.0)) {
          shift = sqrt(b*b + c);
          if (b < 0.0) {
            shift = -shift;
          }
          shift = c/(b + shift);
        }
        double f = (sk + sp)*(sk - sp) + shift;
        double g = sk*ek;
    
        // Chase zeros.
    
        for (j = k; j < p-1; j++) {
          //t = hypot(f,g);
          //Hypotenuse without under/overflow
          t = f*sqrt(1+(g/f)*(g/f));
          double cs = f/t;
          double sn = g/t;
          if (j != k) {
            e.at(j-1) = t;
          }
          f = cs*S.at(j) + sn*e.at(j);
          e.at(j) = cs*e.at(j) - sn*S.at(j);
          g = sn*S.at(j+1);
          S.at(j+1) = cs*S.at(j+1);
          if (wantv) {
            for (i = 0; i < n; i++) {
              t = cs*V.at(i).at(j) + sn*V.at(i).at(j+1);
              V.at(i).at(j+1) = -sn*V.at(i).at(j) + cs*V.at(i).at(j+1);
              V.at(i).at(j) = t;
            }
          }
          //t = hypot(f,g);
          //Hypotenuse without under/overflow
          t = f*sqrt(1+(g/f)*(g/f));
          cs = f/t;
          sn = g/t;
          S.at(j) = t;
          f = cs*e.at(j) + sn*S.at(j+1);
          S.at(j+1) = -sn*e.at(j) + cs*S.at(j+1);
          g = sn*e.at(j+1);
          e.at(j+1) = cs*e.at(j+1);
          if (wantu && (j < m-1)) {
            for (i = 0; i < m; i++) {
              t = cs*U.at(i).at(j) + sn*U.at(i).at(j+1);
              U.at(i).at(j+1) = -sn*U.at(i).at(j) + cs*U.at(i).at(j+1);
              U.at(i).at(j) = t;
            }
          }
        }
        e.at(p-2) = f;
        iter = iter + 1;
      }
      break;
 
      // Convergence.

      case 4: {
 
        // Make the singular values positive.
    
        if (S.at(k) <= 0.0) {
          S.at(k) = (S.at(k) < 0.0 ? -S.at(k) : 0.0);
          if (wantv) {
            for (i = 0; i <= pp; i++) {
              V.at(i).at(k) = -V.at(i).at(k);
            }
          }
        }
    
        // Order the singular values.
    
        while (k < pp) {
          if (S.at(k) >= S.at(k+1)) {
            break;
          }
          t = S.at(k);
          S.at(k) = S.at(k+1);
          S.at(k+1) = t;
          if (wantv && (k < n-1)) {
            for (i = 0; i < n; i++) {
              t = V.at(i).at(k+1);
              V.at(i).at(k+1) = V.at(i).at(k);
              V.at(i).at(k) = t;
            }
          }
          if (wantu && (k < m-1)) {
            for (i = 0; i < m; i++) {
              t = U.at(i).at(k+1);
              U.at(i).at(k+1) = U.at(i).at(k); 
              U.at(i).at(k) = t;
            }
          }
          k++;
        }
        iter = 0;
        p--;
      }
      break;
    } //Switch
  } //While
}
