//Sean M. Law`

#include "LinAlg.hpp"

#include "Misc.hpp"

#include <cmath>
#include <cstdlib>

//modPDB.exe tests/file.pdb -fit tests/file.ref.pdb -outsel :.CA > t.pdb



void LinAlg::SVD (const std::vector< std::vector <double> > &M){
  std::vector< std::vector<double> > A; //A copy of matrix M
  std::vector<double> e;
  std::vector<double> work;
  int m; //Rows
  int n; //Columns
  int last;
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
  
  A=M;
  m=A.size();
  n=A.at(0).size();
  //Check the size of each row
  last=0;
  for (i=0; i<m; i++){
    n=A.at(i).size();
    if (last !=0 && last != n){
      std::cerr << std::endl << "Warning: Matrix has jagged rows" << std::endl;
    }
    last=n;
  }

  nu=std::min(m,n); 
  s.resize(std::min(m+1,n));
  U.resize(m);
  V.resize(n);

  for (i=0; i<m; i++){
    U.at(i).resize(nu);
    for(j=0; j<nu; j++){
      U.at(i).at(j)=0.0;
    }
  }

  for (i=0; i<n; i++){
    V.at(i).resize(n);
    for (j=0; j<n; j++){
      V.at(i).at(j)=0.0;
    }
  }

  e.resize(n);
	for (i=0; i<n; i++){
		e.at(i)=0.0;
	}

  work.resize(m);
	for (i=0; i<m; i++){
		work.at(i)=0.0;
	}

  i=j=k=0;

  // Reduce A to bidiagonal form, storing the diagonal elements
  // in s and the super-diagonal elements in e.
 
  nct = std::min(m-1,n);
  nrt = std::max(0,std::min(n-2,m));
  for (k = 0; k < std::max(nct,nrt); k++) {
    if (k < nct) {
 
      // Compute the transformation for the k-th column and
      // place the k-th diagonal in s.at(k).
      // Compute 2-norm of k-th column without under/overflow.
      s.at(k) = 0;
      for (i = k; i < m; i++) {
        s.at(k) = Misc::hypot(s.at(k),A.at(i).at(k));
      }
      if (s.at(k) != 0.0) {
        if (A.at(k).at(k) < 0.0) {
          s.at(k) = -s.at(k);
        }
        for (i = k; i < m; i++) {
          A.at(i).at(k) /= s.at(k);
        }
        A.at(k).at(k) += 1.0;
      }
      s.at(k) = -s.at(k);
    }
    for (j = k+1; j < n; j++) {
      if ((k < nct) && (s.at(k) != 0.0))  {

        // Apply the transformation.
 
        t = 0;
        for (i = k; i < m; i++) {
          t += A.at(i).at(k)*A.at(i).at(j);
        }
        t = -t/A.at(k).at(k);
        for (i = k; i < m; i++) {
          A.at(i).at(j) += t*A.at(i).at(k);
        }
      }
 
      // Place the k-th row of A into e for the
      // subsequent calculation of the row transformation.
      e.at(j) = A.at(k).at(j);
    }

    if (wantu & (k < nct)) {
 
      // Place the transformation in U for subsequent back
      // multiplication.
 
      for (i = k; i < m; i++) {
        U.at(i).at(k) = A.at(i).at(k);
      }
    }
    if (k < nrt) {
 
      // Compute the k-th row transformation and place the
      // k-th super-diagonal in e.at(k).
      // Compute 2-norm without under/overflow.
      e.at(k) = 0;
      for (i = k+1; i < n; i++) {
        e.at(k) = Misc::hypot(e.at(k),e.at(i));
      }
      if (e.at(k) != 0.0) {
        if (e.at(k+1) < 0.0) {
          e.at(k) = -e.at(k);
        }
        for (i = k+1; i < n; i++) {
           e.at(i) /= e.at(k);
        }
        e.at(k+1) += 1.0;
      }
      e.at(k) = -e.at(k);
      if ((k+1 < m) & (e.at(k) != 0.0)) {

        // Apply the transformation.
        for (i = k+1; i < m; i++) {
          work.at(i) = 0.0;
        }
        for (j = k+1; j < n; j++) {
           for (i = k+1; i < m; i++) {
              work.at(i) += e.at(j)*A.at(i).at(j);
           }
        }
        for (j = k+1; j < n; j++) {
          t = -e.at(j)/e.at(k+1);
          for (i = k+1; i < m; i++) {
            A.at(i).at(j) += t*work.at(i);
          }
        }
      }

      if (wantv) {
 
        // Place the transformation in V for subsequent
        // back multiplication.
 
        for (i = k+1; i < n; i++) {
          V.at(i).at(k) = e.at(i);
        }
      }
    }
  }

  // Set up the final bidiagonal matrix or order p.
 
  p = std::min(n,m+1);
  if (nct < n) {
    s.at(nct) = A.at(nct).at(nct);
  }
  if (m < p) {
    s.at(p-1) = 0.0;
  }
  if (nrt+1 < p) {
    e.at(nrt) = A.at(nrt).at(p-1);
  }
  e.at(p-1) = 0.0;

  // If required, generate U.

  if (wantu) {
    for (j = nct; j < nu; j++) {
      for (i = 0; i < m; i++) {
        U.at(i).at(j) = 0.0;
      }
      U.at(j).at(j) = 1.0;
    }
    for (k = nct-1; k >= 0; k--) {
      if (s.at(k) != 0.0) {
        for (j = k+1; j < nu; j++) {
          t = 0;
          for (i = k; i < m; i++) {
            t += U.at(i).at(k)*U.at(i).at(j);
          }
          t = -t/U.at(k).at(k);
          for (i = k; i < m; i++) {
            U.at(i).at(j) += t*U.at(i).at(k);
          }
        }
        for (i = k; i < m; i++ ) {
          U.at(i).at(k) = -U.at(i).at(k);
        }
        U.at(k).at(k) = 1.0 + U.at(k).at(k);
        for (i = 0; i < k-1; i++) {
          U.at(i).at(k) = 0.0;
        }
      }
      else {
        for (i = 0; i < m; i++) {
          U.at(i).at(k) = 0.0;
        }
        U.at(k).at(k) = 1.0;
      }
    }
  }
 
  // If required, generate V.
  if (wantv) {
    for (k = n-1; k >= 0; k--) {
      if ((k < nrt) & (e.at(k) != 0.0)) {
        for (j = k+1; j < nu; j++) {
          t = 0;
          for (i = k+1; i < n; i++) {
            t += V.at(i).at(k)*V.at(i).at(j);
          }
          t = -t/V.at(k+1).at(k);
          for (i = k+1; i < n; i++) {
            V.at(i).at(j) += t*V.at(i).at(k);
          }
        }
      }
      for (i = 0; i < n; i++) {
        V.at(i).at(k) = 0.0;
      }
      V.at(k).at(k) = 1.0;
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
 
    // kase = 1     if s(p) and e.at(k-1) are negligible and k<p
    // kase = 2     if s(k) is negligible and k<p
    // kase = 3     if e.at(k-1) is negligible, k<p, and
    //              s(k), ..., s(p) are not negligible (qr step).
    // kase = 4     if e(p-1) is negligible (convergence).
 
    for (k = p-2; k >= -1; k--) {
      if (k == -1) {
        break;
      }
      if (abs(e.at(k)) <= eps*(abs(s.at(k)) + abs(s.at(k+1)))) {
        e.at(k) = 0.0;
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
        t = (ks != p ? abs(e.at(ks)) : 0.) + (ks != k+1 ? abs(e.at(ks-1)) : 0.);
        if (abs(s.at(ks)) <= eps*t)  {
          s.at(ks) = 0.0;
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
          t = Misc::hypot(s.at(j),f);
          double cs = s.at(j)/t;
          double sn = f/t;
          s.at(j) = t;
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
          t = Misc::hypot(s.at(j),f);
          double cs = s.at(j)/t;
          double sn = f/t;
          s.at(j) = t;
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
                        abs(s.at(p-1)),abs(s.at(p-2))),abs(e.at(p-2))), 
                        abs(s.at(k))),abs(e.at(k)));
        double sp = s.at(p-1)/scale;
        double spm1 = s.at(p-2)/scale;
        double epm1 = e.at(p-2)/scale;
        double sk = s.at(k)/scale;
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
          double t = Misc::hypot(f,g);
          double cs = f/t;
          double sn = g/t;
          if (j != k) {
            e.at(j-1) = t;
          }
          f = cs*s.at(j) + sn*e.at(j);
          e.at(j) = cs*e.at(j) - sn*s.at(j);
          g = sn*s.at(j+1);
          s.at(j+1) = cs*s.at(j+1);
          if (wantv) {
            for (i = 0; i < n; i++) {
              t = cs*V.at(i).at(j) + sn*V.at(i).at(j+1);
              V.at(i).at(j+1) = -sn*V.at(i).at(j) + cs*V.at(i).at(j+1);
              V.at(i).at(j) = t;
            }
          }
          t = Misc::hypot(f,g);
          cs = f/t;
          sn = g/t;
          s.at(j) = t;
          f = cs*e.at(j) + sn*s.at(j+1);
          s.at(j+1) = -sn*e.at(j) + cs*s.at(j+1);
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
   
        if (s.at(k) <= 0.0) {
          s.at(k) = (s.at(k) < 0.0 ? -s.at(k) : 0.0);
          if (wantv) {
            for (i = 0; i <= pp; i++) {
              V.at(i).at(k) = -V.at(i).at(k);
            }
          }
        }

        // Order the singular values.
   
        while (k < pp) {
          if (s.at(k) >= s.at(k+1)) {
            break;
          }
          t = s.at(k);
          s.at(k) = s.at(k+1);
          s.at(k+1) = t;
          if (wantv && (k < n-1)) {
            for (i = 0; i < n; i++) {
              t = V.at(i).at(k+1); V.at(i).at(k+1) = V.at(i).at(k); V.at(i).at(k) = t;
            }
          }
          if (wantu && (k < m-1)) {
            for (i = 0; i < m; i++) {
              t = U.at(i).at(k+1); U.at(i).at(k+1) = U.at(i).at(k); U.at(i).at(k) = t;
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
