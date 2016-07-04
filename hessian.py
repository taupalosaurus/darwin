

computeEigVal_str = """

  double met[4] = {hess[0], 0.5*(hess[1]+hess[2]), hess[3]};

  int    i;
  double eigVal[2], eigVec[4], nrm,sign,vecNrm,dd,tmp;
  
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


computeEigVal_str2 = """

  double met[4] = {hess[0], 0.5*(hess[1]+hess[2]), hess[3]};

  int    i;
  double eigVal[2], eigVec[4], nrm,sign,vecNrm,dd,tmp;
  
  nrm = fabs(met[0]);  
  for (i=1; i<3; ++i) {
    if ( fabs(met[i]) > nrm )  nrm = fabs(met[i]);
  }    
    
  //--- a null matrix
  if ( nrm < 1e-100 ) {
    eigVal[0] = eigVal[1]= 0.;
    eigVec[0] = 1.; eigVec[1] = 0.;
    eigVec[2] = 0.; eigVec[3] = 1.;
    goto end;
  }
  
  //--- compute eigenvalues
  
  dd = (met[0]-met[2])*(met[0]-met[2]) + 4.*met[1]*met[1];

  //--- Diagonal matrix with two identical eigenvalues
  if ( fabs(dd) < 1.e-24 ) {  // eigVal[0] = eigVal[1] = nrm (or 1 after normalization)
    eigVal[0] = eigVal[1] = met[0];
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
	    eigVal[0] = eigVal[1] = met[0];
        eigVec[0] = 1.; eigVec[1] = 0.;
        eigVec[2] = 0.; eigVec[3] = 1.;
	    goto end;
      }
      else
        exit(13);
    }
  }
  
//  vecNrm     = 1./sqrt(vecNrm);
//  eigVec[2] *= vecNrm;
//  eigVec[3] *= vecNrm;

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

  eigVec[0] =  1;
  eigVec[1] =  0;
  eigVec[2] = 0;
  eigVec[3] =  1;
  eigVal[0] = 1;
  eigVal[1] = 1;

//  eigVal[0] *= nrm;
//  eigVal[1] *= nrm;

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