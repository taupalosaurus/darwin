

computeEigVal_str = """

  double met[4] = {hess[0], 0.5*(hess[1]+hess[2]), hess[3]};

  int    i;
  double eigVal[2], eigVec[4], nrm, sign, vecNrm, dd, tmp;
  
  nrm = fabs(met[0]);
  if ( met[0] >= 0. ) sign =  1.;
  else sign = -1;  
  for (i=1; i<3; ++i) {
    if ( fabs(met[i]) > nrm ) {
      nrm = fabs(met[i]);
      if ( met[i] >= 0. ) sign =  1.;
      else  sign = -1;  
    }  
  }    
    
  //--- a null matrix
  if ( nrm < 1e-100 ) {
    eigVal[0] = eigVal[1]= 0.;
    eigVec[0] = 1.; eigVec[1] = 0.;
    eigVec[2] = 0.; eigVec[3] = 1.;
    goto end;
  }
  
  nrm = sign*nrm;
  
  //--- compute eigenvalues
  for (i=0; i<3; ++i) met[i] = met[i]/nrm;
  
  dd = (met[0]-met[2])*(met[0]-met[2]) + 4.*met[1]*met[1];

  //--- Diagonal matrix with two identical eigenvalues
  if ( fabs(dd) < 1.e-24 ) {  // eigVal[0] = eigVal[1] = nrm (or 1 after normalization)
    eigVal[0] = eigVal[1] = nrm;
    eigVec[0] = 1.; eigVec[1] = 0.;
    eigVec[2] = 0.; eigVec[3] = 1.;
    goto end;
  }
  else if ( dd < 0. ) exit(10);
  
  dd = sqrt(dd);
  eigVal[0] = 0.5*(met[0]+met[2]-dd);
  eigVal[1] = 0.5*(met[0]+met[2]+dd);
  
  //--- compute eigenvectors
  eigVec[2] = -met[1];
  eigVec[3] =  met[0]-eigVal[1];
  vecNrm    =  eigVec[2]*eigVec[2] + eigVec[3]*eigVec[3];

  if ( vecNrm < 1e-30 ) { // => diag matrix => bad sol use the other
    eigVec[2] = eigVal[1]-met[2];
    eigVec[3] = met[1];
    vecNrm    = eigVec[2]*eigVec[2] + eigVec[3]*eigVec[3];
    if ( vecNrm < 1e-30 ) {	//--- in the case we have dd = 0 ie egv1 = egv2 
	  //--- thus M is the Id matrix after normalisation
      if ( fabs(eigVal[0]-1.) < 1.e-12 && fabs(eigVal[1]-1.) < 1.e-12 ) {  
	    eigVal[0] = eigVal[1] = nrm;
        eigVec[0] = 1.; eigVec[1] = 0.;
        eigVec[2] = 0.; eigVec[3] = 1.;
	    goto end;
      }
      else
        exit(13);
    }
  }
  
  vecNrm     = 1./sqrt(vecNrm);
  eigVec[2] *= vecNrm;
  eigVec[3] *= vecNrm;

  if ( fabs(eigVal[0]) < fabs(eigVal[1]) ) {
    eigVec[0] =  eigVec[3];
    eigVec[1] = -eigVec[2];
  }
  else {
    eigVec[0] =  eigVec[2];
    eigVec[1] =  eigVec[3];
    eigVec[2] = -eigVec[1];
    eigVec[3] =  eigVec[0];
    tmp       = eigVal[0];
    eigVal[0] = eigVal[1];
    eigVal[1] = tmp;
  }

  eigVal[0] *= nrm;
  eigVal[1] *= nrm;

  end:
"""



rebuildHessian_str = """

  hess[0] = eigVal[0]*eigVec[0]*eigVec[0] + eigVal[1]*eigVec[2]*eigVec[2];
  hess[1] = eigVal[0]*eigVec[0]*eigVec[1] + eigVal[1]*eigVec[2]*eigVec[3];
  hess[2] = hess[1];
  hess[3] = eigVal[0]*eigVec[1]*eigVec[1] + eigVal[1]*eigVec[3]*eigVec[3];

"""


absValueHessian_str = """

  eigVal[0] = fabs(eigVal[0]);
  eigVal[1] = fabs(eigVal[1]);

"""


truncLowHessian_str = """

  eigVal[0] = fmax(eigVal[0], lmin);
  eigVal[1] = fmax(eigVal[1], lmin);

"""

truncHighHessian_str = """

  eigVal[0] = fmin(eigVal[0], lmax);
  eigVal[1] = fmin(eigVal[1], lmax);

"""


truncRatioHessian_str = """

    double maxLbd = fmax(eigVal[0], eigVal[1]);
    eigVal[0] = fmax(eigVal[0], usa2*maxLbd);
    eigVal[1] = fmax(eigVal[1], usa2*maxLbd);

"""







computeEigVal3_str = """

  int    i, k, cas;
  double eigVal[3], eigVec[9], nrm;
  double us6, us3, a, b, c, d, ap, bp, cp, alpha, beta, lbd1, lbd2, lbd3, delta, eps;
  double px, ppx, x0, x1, xmin, pmin, p[3], tmp; 
  double w1[3], w2[3], w3[3], v1[3], v2[3], v3[3];
  double v1nrm, v2nrm, v3nrm, w1nrm, w2nrm, w3nrm, vecNrm, nrmInv;
  
  
  us6 = 1./6;
  us3 = 1./3;

  double mat[6] = {hess[0], 0.5*(hess[1]+hess[3]), 0.5*(hess[2]+hess[6]), 
                            hess[4],               0.5*(hess[5]+hess[7]),
                                                   hess[8]};
  
  nrm = fabs(mat[0]);
  for (i=1; i<6; ++i)
    nrm = fmax(nrm,fabs(mat[i]));
  
  // check for null matrices
  if ( nrm < 1e-100 ) {
    eigVal[0] = 0; eigVal[1] = 0; eigVal[2] = 0;
    for (i=0; i<9; ++i) eigVec[i] = 0;
    eigVec[0] = 1; eigVec[4] = 1; eigVec[8] = 1;  
    goto end;
  }
  
  // normalize the matrix
  nrmInv = 1. / nrm;
  for (i=0; i<6; ++i) mat[i] *= nrmInv;
  

  a = -mat[0]-mat[3]-mat[5];			// = -Trace(mat)
  b = mat[0]*mat[3] + mat[0]*mat[5] + mat[3]*mat[5]
    - mat[1]*mat[1] - mat[2]*mat[2] - mat[4]*mat[4];
  c = mat[0]*(mat[4]*mat[4]-mat[3]*mat[5])	// = -Det(mat)
    + mat[1]*(mat[1]*mat[5]-mat[2]*mat[4])
    + mat[2]*(mat[2]*mat[3]-mat[1]*mat[4]);
  
  cas = 0;
  
  // P'(X) = ap*X^2 + bp*X + cp
  ap = 3;
  bp = 2*a;
  cp = b;
  
  // First look for double or triple roots  
  delta = bp*bp-4*ap*cp;
  eps   = bp*bp*1e-10;
  
  if ( delta > eps ) {	 
    
    // P' has two different roots:
    // if a root of P' is a root of P, it is a double root for P
    delta  = sqrt(delta);
    eigVal[0] = (-bp+delta)*us6; // first root of P'
    px = (((eigVal[0]+a)*eigVal[0])+b)*eigVal[0] +c;	
    if ( fabs(px) < 1e-15 ) {
      eigVal[1] = eigVal[0];
      eigVal[2] = -a - 2*eigVal[0]; // sum of eig. val. = Trace  (and -a = Trace)
      cas = 2;
      goto endEigVal;
    }
    eigVal[0] = (-bp-delta)*us6; // other root of P' 
    px = (((eigVal[0]+a)*eigVal[0])+b)*eigVal[0] +c;
    if ( fabs(px) < 1e-15 ) {
      eigVal[1] = eigVal[0];
      eigVal[2] = -a - 2*eigVal[0];
      cas = 2;
      goto endEigVal;
    }
    // else P has 3 single roots, see later
  }
  else if ( fabs(delta) <= eps ) { // delta = 0 
    // P' has a double root => P has a triple root    
    eigVal[0] = eigVal[1] = eigVal[2] = -bp*us6;
    cas = 3;
    goto endEigVal;
  }
  else {
    exit(11);   // P' cannot have no real root
  }
  
  // If P has 3 single roots, find the middle one with a Newton algorithm
  // the inflection point (ie P"(x)=0) is used as a starting point
  // TODO is this faster than using Cardano's formulas ?
  
  x0  = -a*us3;               // root of P''
  ppx = b-a*a*us3;            // P'(x0) 
  px = (2*a*a*us6+b)*a*us3+c;	// P(x0)
  
  xmin = x0;
  pmin = fabs(px);
  
  for (i=0; i<100; ++i) {
    
    x1  = x0-px/ppx;
    px = (((x1+a)*x1)+b)*x1+c;
    if ( fabs(px) < 1e-18 ) {
      eigVal[1] = x1;
      break;
    }
    if ( fabs(px) < pmin ) {
      xmin = x1;
      pmin = fabs(px);
    }
    ppx = (ap*x1+bp)*x1+cp; 
    if ( fabs((x1-x0)/x1) < 1e-7 ) {  // the algorithm has (almost) converged
      eigVal[1] = xmin;
      if ( pmin > 1e-12 ){ // check that P(xmin) really close to 0
	      printf(\"ERRORKERNEL pmin: %1.3e, dist: %1.2e i: %d    mat: %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\\n\", pmin, fabs((x1-x0)/x1), i, mat[0], mat[1], mat[2], mat[3], mat[4], mat[5]);
        printf(\"ERRORKERNEL  eigVals so far: %1.8e  %1.8e %1.8e\\n\", eigVal[0], eigVal[1], eigVal[2]);
        exit(22);
      }
      break;
    } 
    x0 = x1;
  }
  
  if ( i == 100 ) {
    eigVal[1] = xmin;
    if ( pmin > 1e-12 )
      exit(13);
  }
  
  
  // We have now P(X) = (X-eigVal[1])*(X^2 + alpha*X + beta)
  // so we just have to find the roots of quadratic polynomial
  alpha = a + eigVal[1];
  beta = b + eigVal[1]*alpha;
  
  delta = alpha*alpha-4*beta;
  eps   = alpha*alpha*1e-10; 
  
  if ( fabs(delta) < eps ) { // double root (lamda[0]=eigVal[1])
    eigVal[2] = eigVal[1];
    eigVal[0] = eigVal[1] = -0.5*alpha;
    cas = 2;
    goto endEigVal;
  }
  else if ( delta < 0. ) {
    exit(14);
  }
  
  delta  = sqrt(delta);
  eigVal[2] =  0.5*(delta-alpha); 
  eigVal[0] = -0.5*(delta+alpha); 
  
  //--- check another time if a double/triple roots is obtained at 1e-5
  //--- very important for the conditionning when looking for the eigenvectors
  for (i=0; i<3; ++i)
    p[i] = fabs( (((eigVal[i]+a)*eigVal[i] )+b)*eigVal[i]+c);
  lbd1 = fabs(eigVal[0]);
  lbd2 = fabs(eigVal[1]);
  lbd3 = fabs(eigVal[2]);
  
  if ( fabs(eigVal[0]-eigVal[1]) < 1e-5*(lbd1+lbd2)*0.5 ) {
    if ( p[0] < p[1] ) {
      eigVal[1] = eigVal[0];
      p[1] = p[0];
    }
    else {
      eigVal[0] = eigVal[1];
      p[0] = p[1];
    }
      
    // check for triple root
    if ( fabs(eigVal[1]-eigVal[2]) < 1e-5*(lbd2+lbd3)*0.5 ) {
      if ( p[1] < p[2] ) eigVal[2] = eigVal[1];
      else eigVal[1] = eigVal[2];
      cas = 3;
      goto endEigVal;
    }
    else 
      cas = 2;
      goto endEigVal;
  }
  
  else if ( fabs(eigVal[1]-eigVal[2]) < 1e-5*(lbd2+lbd3)*0.5 ) {
    if ( p[2] < p[1] ) {
      eigVal[1] = eigVal[2];
      p[1] = p[2];
    }
    else {
      eigVal[2] = eigVal[1];
      p[2] = p[1];
    }
    tmp = eigVal[0];	
    eigVal[0] = eigVal[2];
    eigVal[2] = tmp;    
    cas = 2;
    goto endEigVal;
  }
  
  else if ( fabs(eigVal[0]-eigVal[2]) < 1e-5*(lbd1+lbd3)*0.5 ) {
    if ( p[2] < p[0] ) {
      eigVal[0] = eigVal[2];
      p[0] = p[2];
    }
    else {
      eigVal[2] = eigVal[0];
      p[2] = p[0];
    }
    tmp = eigVal[1];
    eigVal[1] = eigVal[2];
    eigVal[2] = tmp;    
    cas = 2;
    goto endEigVal;
  }
  
  cas = 1;
  endEigVal:
  
  switch (cas) {
  
  // case with hree single eigenValues
  case 1:
    // vk= wi/\wj !=0  with  wi,wj line of W = mat-eigVal[i]*Id
    for (k=0; k<2; ++k) {
      w1[0] = mat[0]-eigVal[k]; w1[1] = mat[1];           w1[2] = mat[2];
      w2[0] = mat[1];           w2[1] = mat[3]-eigVal[k]; w2[2] = mat[4];
      w3[0] = mat[2];           w3[1] = mat[4];           w3[2] = mat[5]-eigVal[k];

      // cross products
      v1[0] =  w1[1]*w3[2] - w1[2]*w3[1]; v1[1] = -w1[0]*w3[2] + w1[2]*w3[0]; v1[2] =  w1[0]*w3[1] - w1[1]*w3[0];
      v1nrm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
      v2[0] =  w1[1]*w2[2] - w1[2]*w2[1]; v2[1] = -w1[0]*w2[2] + w1[2]*w2[0]; v2[2] =  w1[0]*w2[1] - w1[1]*w2[0];
      v2nrm = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2];
      v3[0] =  w2[1]*w3[2] - w2[2]*w3[1]; v3[1] = -w2[0]*w3[2] + w2[2]*w3[0]; v3[2] =  w2[0]*w3[1] - w2[1]*w3[0];
      v3nrm = v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2];

      // take the vector with the highest norm
      if ( v1nrm >= v2nrm && v1nrm >= v3nrm) {
       vecNrm = 1./sqrt(v1nrm);
        eigVec[3*k] = v1[0]*vecNrm; eigVec[3*k+1] = v1[1]*vecNrm; eigVec[3*k+2] = v1[2]*vecNrm;
      }
      else if ( v2nrm >= v1nrm && v2nrm >= v3nrm) {
       vecNrm = 1./sqrt(v2nrm);
        eigVec[3*k] = v2[0]*vecNrm; eigVec[3*k+1] = v2[1]*vecNrm; eigVec[3*k+2] = v2[2]*vecNrm;
      }
      else {
       vecNrm = 1./sqrt(v3nrm);
        eigVec[3*k] = v3[0]*vecNrm; eigVec[3*k+1] = v3[1]*vecNrm; eigVec[3*k+2] = v3[2]*vecNrm;
      }
    }

    // The last eigenvector is simply orthogonal to both others : v3=v1/\v2
    eigVec[6]  = eigVec[1]*eigVec[5] - eigVec[2]*eigVec[4];
    eigVec[7]  = eigVec[2]*eigVec[3] - eigVec[0]*eigVec[5];
    eigVec[8]  = eigVec[0]*eigVec[4] - eigVec[1]*eigVec[3];
   vecNrm = eigVec[6]*eigVec[6]+eigVec[7]*eigVec[7]+eigVec[8]*eigVec[8];
   vecNrm = 1./sqrt(vecNrm);
    eigVec[6] *=vecNrm; 
    eigVec[7] *=vecNrm; 
    eigVec[8] *=vecNrm; 
    
    break;

  // case with a double eigVal (eigVal[0] = eigVal[1]) and a single one (eigVal[2])  
  case 2:
    //  v3
    w1[0] = mat[0]-eigVal[2]; w1[1] = mat[1];           w1[2] = mat[2];
    w2[0] = mat[1];           w2[1] = mat[3]-eigVal[2]; w2[2] = mat[4];
    w3[0] = mat[2];           w3[1] = mat[4];           w3[2] = mat[5]-eigVal[2];

    // cross product of the lines
    v1[0] =  w1[1]*w3[2] - w1[2]*w3[1]; v1[1] = -w1[0]*w3[2] + w1[2]*w3[0]; v1[2] =  w1[0]*w3[1] - w1[1]*w3[0];
    v1nrm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
    v2[0] =  w1[1]*w2[2] - w1[2]*w2[1]; v2[1] = -w1[0]*w2[2] + w1[2]*w2[0]; v2[2] =  w1[0]*w2[1] - w1[1]*w2[0];
    v2nrm = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2];
    v3[0] =  w2[1]*w3[2] - w2[2]*w3[1]; v3[1] = -w2[0]*w3[2] + w2[2]*w3[0]; v3[2] =  w2[0]*w3[1] - w2[1]*w3[0];
    v3nrm = v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2];

    // take the vector with the highest norm
    if ( v1nrm >= v2nrm && v1nrm >= v3nrm) {
     vecNrm = 1./sqrt(v1nrm);
      eigVec[6] = v1[0]*vecNrm; eigVec[7] = v1[1]*vecNrm; eigVec[8] = v1[2]*vecNrm;
    }
    else if ( v2nrm >= v1nrm && v2nrm >= v3nrm) {
     vecNrm = 1./sqrt(v2nrm);
      eigVec[6] = v2[0]*vecNrm; eigVec[7] = v2[1]*vecNrm; eigVec[8] = v2[2]*vecNrm;
    }
    else {
     vecNrm = 1./sqrt(v3nrm);
      eigVec[6] = v3[0]*vecNrm; eigVec[7] = v3[1]*vecNrm; eigVec[8] = v3[2]*vecNrm;
    }    
    
    // v1
    w1nrm = w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2];
    w2nrm = w2[0]*w2[0]+w2[1]*w2[1]+w2[2]*w2[2];
    w3nrm = w3[0]*w3[0]+w3[1]*w3[1]+w3[2]*w3[2];
    if ( w1nrm >=  w2nrm && w1nrm > w3nrm ) {
     vecNrm = 1./sqrt(w1nrm);
      eigVec[0] = w1[0]*vecNrm; eigVec[1] = w1[1]*vecNrm; eigVec[2] = w1[2]*vecNrm;
    }
    else if ( w2nrm >=  w1nrm && w2nrm > w3nrm ) {
     vecNrm = 1./sqrt(w2nrm);
      eigVec[0] = w2[0]*vecNrm; eigVec[1] = w2[1]*vecNrm; eigVec[2] = w2[2]*vecNrm;
    }
    else {
     vecNrm = 1./sqrt(w3nrm);
      eigVec[0] = w3[0]*vecNrm; eigVec[1] = w3[1]*vecNrm; eigVec[2] = w3[2]*vecNrm;
    }
        
    // The last eigenvector is simply orthogonal to both others: v2=v3/\v1
    eigVec[3]  = eigVec[7]*eigVec[2] - eigVec[8]*eigVec[1];
    eigVec[4]  = eigVec[8]*eigVec[0] - eigVec[6]*eigVec[2];
    eigVec[5]  = eigVec[6]*eigVec[1] - eigVec[7]*eigVec[0];
   vecNrm = eigVec[3]*eigVec[3]+eigVec[4]*eigVec[4]+eigVec[5]*eigVec[5];
   vecNrm = 1./sqrt(vecNrm);
    eigVec[3] *=vecNrm; eigVec[4] *=vecNrm; eigVec[5] *=vecNrm; 
    break;

  // triple eigenvalue => isotropic metric = lamba*Id  
  case 3:
    for (i=0; i<9; ++i) 
      eigVec[i] = 0;
    eigVec[0] = 1; eigVec[4] = 1; eigVec[8] = 1;  
    break;  
  
  default:
    exit(15);
  }
  
  eigVal[0] *= nrm;
  eigVal[1] *= nrm;
  eigVal[2] *= nrm;
  
  end:

"""


rebuildHessian3_str = """
  
  hess[0] = eigVal[0]*eigVec[0]*eigVec[0] + eigVal[1]*eigVec[3]*eigVec[3] + eigVal[2]*eigVec[6]*eigVec[6];
  hess[1] = eigVal[0]*eigVec[0]*eigVec[1] + eigVal[1]*eigVec[3]*eigVec[4] + eigVal[2]*eigVec[6]*eigVec[7];
  hess[2] = eigVal[0]*eigVec[0]*eigVec[2] + eigVal[1]*eigVec[3]*eigVec[5] + eigVal[2]*eigVec[6]*eigVec[8];
  hess[3] = hess[1];
  hess[4] = eigVal[0]*eigVec[1]*eigVec[1] + eigVal[1]*eigVec[4]*eigVec[4] + eigVal[2]*eigVec[7]*eigVec[7];
  hess[5] = eigVal[0]*eigVec[1]*eigVec[2] + eigVal[1]*eigVec[4]*eigVec[5] + eigVal[2]*eigVec[7]*eigVec[8];
  hess[6] = hess[2];
  hess[7] = hess[5];
  hess[8] = eigVal[0]*eigVec[2]*eigVec[2] + eigVal[1]*eigVec[5]*eigVec[5] + eigVal[2]*eigVec[8]*eigVec[8];
  

"""

absValueHessian3_str = """

  eigVal[0] = fabs(eigVal[0]);
  eigVal[1] = fabs(eigVal[1]);
  eigVal[2] = fabs(eigVal[2]);

"""


truncLowHessian3_str = """

  eigVal[0] = fmax(eigVal[0], lmin);
  eigVal[1] = fmax(eigVal[1], lmin);
  eigVal[2] = fmax(eigVal[2], lmin);

"""

truncHighHessian3_str = """

  eigVal[0] = fmin(eigVal[0], lmax);
  eigVal[1] = fmin(eigVal[1], lmax);
  eigVal[2] = fmin(eigVal[2], lmax);

"""


truncRatioHessian3_str = """

    double maxLbd = fmax(eigVal[0], fmax(eigVal[1], eigVal[2]));
    eigVal[0] = fmax(eigVal[0], usa2*maxLbd);
    eigVal[1] = fmax(eigVal[1], usa2*maxLbd);
    eigVal[2] = fmax(eigVal[2], usa2*maxLbd);

"""


