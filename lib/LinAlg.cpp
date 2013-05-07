//Sean M. Law`

#include "LinAlg.hpp"
//A = [[168422, 52018.2, 879.072],[47630.9,  154571, 23175.8],[411.489, 17325.3, 94222.3]]
//modPDB.exe tests/file.pdb -fit tests/file.ref.pdb -outsel :.CA > t.pdb
void LinAlg::SVD (const std::vector< std::vector <double> > &M){
  std::vector< std::vector<double> > A; //A copy of matrix M
  std::vector< std::vector<double> > U;
  std::vector< std::vector<double> > V;
  std::vector<double> s;
  std::vector<double> e;
  std::vector<double> work;
  int m; //Rows
  int n; //Columns
  int nu;
  bool wantu;
  bool wantv;
  double t;

  int i;
  int j;
  int k;

  int nct;
  int nrt;
  int p;

  int pp;
  int iter;
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
  for (i=0; i< A.size(); i++){
    A.at(i).resize(n);
    for (j=0; j< A.at(i).size(); j++){
      A.at(i).at(j)=M.at(i).at(j);
    }
  }

  nu=std::min(m,n); 
  s.resize(std::min(m+1,n));
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
  i=j=k=0;

  // Reduce A to bidiagonal form, storing the diagonal elements
  // in s and the super-diagonal elements in e.
 
  nct = std::min(m-1,n);
  nrt = std::max(0,std::min(n-2,m));
  for (k = 0; k < std::max(nct,nrt); k++) {
    if (k < nct) {
 
      // Compute the transformation for the k-th column and
      // place the k-th diagonal in s[k].
      // Compute 2-norm of k-th column without under/overflow.
      s[k] = 0;
      for (i = k; i < m; i++) {
        s[k] = Misc::hypot(s[k],A[i][k]);
      }
      if (s[k] != 0.0) {
        if (A[k][k] < 0.0) {
          s[k] = -s[k];
        }
        for (i = k; i < m; i++) {
          A[i][k] /= s[k];
        }
        A[k][k] += 1.0;
      }
      s[k] = -s[k];
    }
    for (j = k+1; j < n; j++) {
      if ((k < nct) && (s[k] != 0.0))  {

        // Apply the transformation.
 
        t = 0;
        for (i = k; i < m; i++) {
          t += A[i][k]*A[i][j];
        }
        t = -t/A[k][k];
        for (i = k; i < m; i++) {
          A[i][j] += t*A[i][k];
        }
      }
 
      // Place the k-th row of A into e for the
      // subsequent calculation of the row transformation.
 
      e[j] = A[k][j];
    }
    if (wantu & (k < nct)) {
 
      // Place the transformation in U for subsequent back
      // multiplication.
 
      for (i = k; i < m; i++) {
        U[i][k] = A[i][k];
      }
    }
    if (k < nrt) {
 
      // Compute the k-th row transformation and place the
      // k-th super-diagonal in e[k].
      // Compute 2-norm without under/overflow.
      e[k] = 0;
      for (i = k+1; i < n; i++) {
        e[k] = Misc::hypot(e[k],e[i]);
      }
      if (e[k] != 0.0) {
        if (e[k+1] < 0.0) {
          e[k] = -e[k];
        }
        for (i = k+1; i < n; i++) {
           e[i] /= e[k];
        }
        e[k+1] += 1.0;
      }
      e[k] = -e[k];
      if ((k+1 < m) & (e[k] != 0.0)) {

        // Apply the transformation.
 
        for (i = k+1; i < m; i++) {
          work[i] = 0.0;
        }
        for (j = k+1; j < n; j++) {
           for (i = k+1; i < m; i++) {
              work[i] += e[j]*A[i][j];
           }
        }
        for (j = k+1; j < n; j++) {
          t = -e[j]/e[k+1];
          for (i = k+1; i < m; i++) {
            A[i][j] += t*work[i];
          }
        }
      }
      if (wantv) {
 
        // Place the transformation in V for subsequent
        // back multiplication.
 
        for (i = k+1; i < n; i++) {
          V[i][k] = e[i];
        }
      }
    }
  }
 
  // Set up the final bidiagonal matrix or order p.
 
  p = std::min(n,m+1);
  if (nct < n) {
    s[nct] = A[nct][nct];
  }
  if (m < p) {
    s[p-1] = 0.0;
  }
  if (nrt+1 < p) {
    e[nrt] = A[nrt][p-1];
  }
  e[p-1] = 0.0;

  // If required, generate U.

  if (wantu) {
    for (j = nct; j < nu; j++) {
      for (i = 0; i < m; i++) {
        U[i][j] = 0.0;
      }
      U[j][j] = 1.0;
    }
    for (k = nct-1; k >= 0; k--) {
      if (s[k] != 0.0) {
        for (j = k+1; j < nu; j++) {
          t = 0;
          for (i = k; i < m; i++) {
            t += U[i][k]*U[i][j];
          }
          t = -t/U[k][k];
          for (i = k; i < m; i++) {
            U[i][j] += t*U[i][k];
          }
        }
        for (i = k; i < m; i++ ) {
          U[i][k] = -U[i][k];
        }
        U[k][k] = 1.0 + U[k][k];
        for (i = 0; i < k-1; i++) {
          U[i][k] = 0.0;
        }
      }
      else {
        for (i = 0; i < m; i++) {
          U[i][k] = 0.0;
        }
        U[k][k] = 1.0;
      }
    }
  }
 
  // If required, generate V.
  if (wantv) {
    for (k = n-1; k >= 0; k--) {
      if ((k < nrt) & (e[k] != 0.0)) {
        for (j = k+1; j < nu; j++) {
          t = 0;
          for (i = k+1; i < n; i++) {
            t += V[i][k]*V[i][j];
          }
          t = -t/V[k+1][k];
          for (i = k+1; i < n; i++) {
            V[i][j] += t*V[i][k];
          }
        }
      }
      for (i = 0; i < n; i++) {
        V[i][k] = 0.0;
      }
      V[k][k] = 1.0;
    }
  }
 
  // Main iteration loop for the singular values.
 
  pp = p-1;
  iter = 0;
  eps = pow(2.0,-52.0);
  while (p > 0) {
    int k=0;
    int kase=0;
 
    // Here is where a test for too many iterations would go.
 
    // This section of the program inspects for
    // negligible elements in the s and e arrays.  On
    // completion the variables kase and k are set as follows.
 
    // kase = 1     if s(p) and e[k-1] are negligible and k<p
    // kase = 2     if s(k) is negligible and k<p
    // kase = 3     if e[k-1] is negligible, k<p, and
    //              s(k), ..., s(p) are not negligible (qr step).
    // kase = 4     if e(p-1) is negligible (convergence).
 
    for (k = p-2; k >= -1; k--) {
      if (k == -1) {
        break;
      }
      if (abs(e[k]) <= eps*(abs(s[k]) + abs(s[k+1]))) {
        e[k] = 0.0;
        break;
      }
    }
    if (k == p-2) {
      kase = 4;
    } 
    else {
      int ks;
      for (ks = p-1; ks >= k; ks--) {
        if (ks == k) {
          break;
        }
        t = (ks != p ? abs(e[ks]) : 0.) + (ks != k+1 ? abs(e[ks-1]) : 0.);
        if (abs(s[ks]) <= eps*t)  {
          s[ks] = 0.0;
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
        double f = e[p-2];
        e[p-2] = 0.0;
        for (j = p-2; j >= k; j--) {
          t = Misc::hypot(s[j],f);
          double cs = s[j]/t;
          double sn = f/t;
          s[j] = t;
          if (j != k) {
            f = -sn*e[j-1];
            e[j-1] = cs*e[j-1];
          }
          if (wantv) {
            for (i = 0; i < n; i++) {
              t = cs*V[i][j] + sn*V[i][p-1];
              V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
              V[i][j] = t;
            }
          }
        }
      }
      break;
 
      // Split at negligible s(k).
      
      case 2: {
        double f = e[k-1];
        e[k-1] = 0.0;
        for (j = k; j < p; j++) {
          t = Misc::hypot(s[j],f);
          double cs = s[j]/t;
          double sn = f/t;
          s[j] = t;
          f = -sn*e[j];
          e[j] = cs*e[j];
          if (wantu) {
            for (i = 0; i < m; i++) {
              t = cs*U[i][j] + sn*U[i][k-1];
              U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
              U[i][j] = t;
            }
          }
        }
      }
      break;
 
      // Perform one qr step.
 
      case 3: {
 
        // Calculate the shift.
    
        double scale = std::max(std::max(std::max(std::max(
                        abs(s[p-1]),abs(s[p-2])),abs(e[p-2])), 
                        abs(s[k])),abs(e[k]));
        double sp = s[p-1]/scale;
        double spm1 = s[p-2]/scale;
        double epm1 = e[p-2]/scale;
        double sk = s[k]/scale;
        double ek = e[k]/scale;
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
          double t = Misc::hypot(f,g);
          double cs = f/t;
          double sn = g/t;
          if (j != k) {
            e[j-1] = t;
          }
          f = cs*s[j] + sn*e[j];
          e[j] = cs*e[j] - sn*s[j];
          g = sn*s[j+1];
          s[j+1] = cs*s[j+1];
          if (wantv) {
            for (i = 0; i < n; i++) {
              t = cs*V[i][j] + sn*V[i][j+1];
              V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
              V[i][j] = t;
            }
          }
          t = Misc::hypot(f,g);
          cs = f/t;
          sn = g/t;
          s[j] = t;
          f = cs*e[j] + sn*s[j+1];
          s[j+1] = -sn*e[j] + cs*s[j+1];
          g = sn*e[j+1];
          e[j+1] = cs*e[j+1];
          if (wantu && (j < m-1)) {
            for (i = 0; i < m; i++) {
              t = cs*U[i][j] + sn*U[i][j+1];
              U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
              U[i][j] = t;
            }
          }
        }
        e[p-2] = f;
        iter = iter + 1;
      }
      break;

      // Convergence.
      
      case 4: {
        // Make the singular values positive.
   
        if (s[k] <= 0.0) {
          s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
          if (wantv) {
            for (i = 0; i <= pp; i++) {
              V[i][k] = -V[i][k];
            }
          }
        }

        // Order the singular values.
   
        while (k < pp) {
          if (s[k] >= s[k+1]) {
            break;
          }
          t = s[k];
          s[k] = s[k+1];
          s[k+1] = t;
          if (wantv && (k < n-1)) {
            for (i = 0; i < n; i++) {
              t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
            }
          }
          if (wantu && (k < m-1)) {
            for (i = 0; i < m; i++) {
              t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
            }
          }
          k++;
        }
        iter = 0;
        p--;
      }
      break;
    }
  }

}
