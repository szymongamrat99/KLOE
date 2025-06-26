#include <cmath>
#include <cstdlib>
#include <ctime>
#include <stdio.h>

#include "reconstructor.h"

#define ITERATORS int _i[6]

#define ITER6CLU( max ) for(_i[0]=0;_i[0]<max-5;_i[0]++){ \
  for(_i[1]=_i[0]+1;_i[1]<max-4;_i[1]++){			\
  for(_i[2]=_i[1]+1;_i[2]<max-3;_i[2]++){			\
  for(_i[3]=_i[2]+1;_i[3]<max-2;_i[3]++){			\
  for(_i[4]=_i[3]+1;_i[4]<max-1;_i[4]++){			\
  for(_i[5]=_i[4]+1;_i[5]<max;_i[5]++)			

#define ENDITER6CLU }}}}}


int combs[15][4] = { {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 3, 6}, 
		     {1, 2, 4, 5}, {1, 2, 4, 6}, {1, 2, 5, 6}, 
		     {1, 3, 4, 5}, {1, 3, 4, 6}, {1, 3, 5, 6}, 
		     {1, 4, 5, 6}, {2, 3, 4, 5}, {2, 3, 4, 6}, 
		     {2, 3, 5, 6}, {2, 4, 5, 6}, {3, 4, 5, 6}};

int leftout[15][2] = { {5,6}, {4,6}, {4,5},
		       {3,6}, {3,5}, {3,4},
		       {2,6}, {2,5}, {2,4},
		       {2,3}, {1,6}, {1,5},
		       {1,4}, {1,3}, {1,2}};

int combs6[15][6] = 
  {{1, 2, 3, 4, 5, 6}, {1, 2, 3, 5, 4, 6}, {1, 2, 3, 6, 4, 5},
   {1, 2, 4, 5, 3, 6}, {1, 2, 4, 6, 3, 5}, {1, 2, 5, 6, 3, 4},
   {1, 3, 4, 5, 2, 6}, {1, 3, 4, 6, 2, 5}, {1, 3, 5, 6, 2, 4},
   {1, 4, 5, 6, 2, 3}, {2, 3, 4, 5, 1, 6}, {2, 3, 4, 6, 1, 5},
   {2, 3, 5, 6, 1, 4}, {2, 4, 5, 6, 1, 3}, {3, 4, 5, 6, 1, 2}};


/********************  Some globals  ********************************/
const double mPi0 = 134.9766; // Mev/c2
const double mK0  = 497.614;  // Mev/c2 

#define ECLMIN 20. // minimal single-cluster energy

/********** checking whether 6-cluster sets are indetical ***********/
bool AreIdentical(int * g6mc, int * set, int * tab){
  for(int i=0;i<6;i++){
    if( g6mc[i] != tab[ set[i] ] ){
      return false;
    }
  }
  return true;
}



/*********  One function to be provided to AIX code  ****************/
extern "C" {

void reconstruct3pi0_(int * ncl,
		      int * ncll,
		      float * Xcl,
		      float * Ycl,
		      float * Zcl,
		      float * Tcl,
		      float * Enecl,
		      int * g6mcok,
		      int * g6mc,
		      // output
		      int * g6taken,
		      float * Knerec,
		      float * KlTime,
		      int * nSetsPassed,
		      int * trueSetFound
		      );

}


void reconstruct3pi0_(int * ncl,
		      int * ncll,
		      float * Xcl,
		      float * Ycl,
		      float * Zcl,
		      float * Tcl,
		      float * Enecl,
		      int * g6mcok,
		      int * g6mc,
		      // output
		      int * g6taken,
		      float * Knerec,
		      float * KlTime,
		      int * nSetsPassed,
		      int * trueSetFound
		      ){
  

  // reset
  int bad = 0;
  
  *nSetsPassed = 0;

  Reconstructor * Rec = new Reconstructor;
  
  /************************ Cluster selection *************************/
  int ngood = 0;
  int bestComb = -1;
  int bestSet[6];
  int combNo = 0; // number of currently considered combination
  
  double mK0Min = 1.0e10;
  double residualMin = 1.0e10;

  double spreadSol[2][4];
  double meanSol[2][4];
  float residual[15][2];   


  ITERATORS;
  ITER6CLU( *ncl ){

    combNo++;
    
    int cc;
    bool clusterEnergyTooLow = false;
    for(int i=0;i<6;i++){
      cc = _i[i];
      cc = ncll[ cc ] - 1;
      // check single-cluster energy criterion
      if( Enecl[cc] < ECLMIN ){
	clusterEnergyTooLow = true;
      }
      // set values into the reconstructor
      Rec->SetClu(i, Xcl[cc], Ycl[cc], Zcl[cc], Tcl[cc], Enecl[cc]);
    }
    
    
    // 1st criterion - single-cluster energies


    if( clusterEnergyTooLow ){
      continue;
    }    


    // 2nd criterion - Energy of 6 clusters
    double totEne = Rec->TotalEnergy();
    if( totEne <  350. || totEne > 650. ){
      continue;
    }    
    
    bad = 0;
    Solution zeroSolution[15];
    double Etemp[2];
    Etemp[0] = 0;
    Etemp[1] = 0;
    for(int i=0;i<2;i++){
      for(int j=0;j<4;j++){
	meanSol[i][j] = 0;
      }
    }     
    double minResidual = 1e10;
    int minComb = -1;
    int ncused[2];
    ncused[0] = 0;
    ncused[1] = 0;
    
    /*************** Iterate over combinations ****************/
    double totalResidual = 0;
    double Cene[15];
    double ng = 0;
      
    for(int ic = 0; ic < 15; ic++){
      zeroSolution[ic] = Rec->MySolve( combs[ic] );
      if( !zeroSolution[ic].error[0] )ng++;
      Cene[ic] = Rec->CombEnergy( combs[ic] );
      // add solution to average
      for(int i=0;i<2;i++){ // over solutions
	if( zeroSolution[ic].error[i] ) continue;
	// only if solution without errors
	for(int j=0;j<4;j++){ // over coordiates
	  meanSol[i][j] +=  zeroSolution[ic].sol[i][j] * Cene[ic];
	}
	Etemp[i] += Cene[ic];
	ncused[i]++;
      } 
      // Calculate residual errors
      residual[ic][0] = 0;
      residual[ic][1] = 0;
	
      for(int i=0; i<2;i++){
	residual[ic][1] += Rec->ResidualErr( leftout[ic][i]-1, zeroSolution[ic].sol[0] );
      }
	
      if( residual[ic][1] == residual[ic][1] ){ // only add to total if not NaN
	totalResidual += residual[ic][1];
	
	// select solution with minimal residual
	if( residual[ic][1] < minResidual ){
	  minComb = ic;
	  minResidual = residual[ic][1];
	}
	
      }
      
    }
    
    if( minComb < 0 ){
      bad = 4;
    }

    /*************** END iterating over combinations ***************/
    
    for(int i=0;i<2;i++){
      for(int j=0;j<4;j++){
	meanSol[i][j] /= Etemp[i];
      }
    }
    // calculate positions spread
    for(int i=0;i<2;i++){
      for(int j=0;j<4;j++){
	spreadSol[i][j] = 0.;
	for(int ic = 0; ic < 15; ic++){
	  if( zeroSolution[ic].error[i] == false){
	    spreadSol[i][j] += Cene[ic]*pow( meanSol[i][j]-zeroSolution[ic].sol[i][j] , 2.);
	  }
	}
	spreadSol[i][j] /= Etemp[i];
	spreadSol[i][j] = sqrt( spreadSol[i][j] );
      }
    }

    
    // get invariant masses
    double errMin = 1.0e10;
    double m2g[3];
    double m6g;
    int icMin = 0;
    
    for(int ic=0;ic<15;ic++){ // over combinations
      m6g = Rec->GetInvMasses( //iS.sol[0]
			      meanSol[0], combs6[ic], m2g );
	
      if( fabs( m2g[0]-mPi0 ) +
	  fabs( m2g[1]-mPi0 ) +
	  fabs( m2g[2]-mPi0 ) < errMin ){
	errMin = ( fabs( m2g[0]-mPi0 ) +
		   fabs( m2g[1]-mPi0 ) +
		   fabs( m2g[2]-mPi0 ) );
	icMin = ic;
      } 
    }
    
    // for best pairing only
    m6g = Rec->GetInvMasses( // iS.sol[0]
			      meanSol[0], combs6[icMin], m2g);

    // 3rd criterion - residual of equations left out (total of all combinations)
    if( totalResidual > 2000. || totalResidual==0 ){
      continue;
    }
    
    // 4th criterion - solutions spread around mean - upper side
    if( spreadSol[0][0] > 100. ||
	spreadSol[0][1] > 100. ||
	spreadSol[0][2] > 100. ||
	spreadSol[0][3] > 4.0 ){
      //      continue;
    }
    
    // 5th criterion - solutions spread around mean - lower side
    if( spreadSol[0][0] < 0.01 ||
	spreadSol[0][1] < 0.01 ||
	spreadSol[0][2] < 0.01 ||
	spreadSol[0][3] < 0.002 ){
      //      continue;
    }
    
    // 6th criterion - reasonability of mean solution
    if( meanSol[0][3] < -10 || meanSol[0][3] > 60.){
      //      continue;
    }else if( pow(meanSol[0][0],2.)+pow(meanSol[0][1],2.) > 190.*190. || fabs( meanSol[0][2] )>170. ){
      //      continue;
    }
    
    if(ngood==0){
      ngood = 1;
      for(int k=0;k<6;k++){
	bestSet[k] = _i[k];
      }
    }
      
    // 7th criterion - 6g invariant mass
    if( m6g < 300. ){
      continue;
    }

    // passed all criteria - count a good cluster set
    *nSetsPassed = *nSetsPassed + 1;
    // fill K momentum in case it was the only good set
    Rec->GetKmomentum( zeroSolution[minComb].sol[0] , Knerec);
    

    
    // find cluster set with smallest deviation from mK0
    //    if( fabs(m6g-mK0) < mK0Min ){
    //      mK0Min = fabs(m6g-mK0) ;
    if( totalResidual < residualMin ){
      residualMin = totalResidual;
      bestComb = combNo;
      // store this set
      for(int k=0;k<6;k++){
	bestSet[k] = _i[k];
      }
      // update best Knerec
      Rec->GetKmomentum( zeroSolution[minComb].sol[0] , Knerec);
      for(int i=0;i<3;i++){
	Knerec[6+i] = zeroSolution[minComb].sol[0][i];
      }
      *KlTime = zeroSolution[minComb].sol[0][3];
    }
    
    ENDITER6CLU}
  /********************** end cluster selection ***********************/
  
  // store best cluster set for output
  for(int k=0;k<6;k++){
    g6taken[k] = bestSet[k];
  }
  
  if( *nSetsPassed > 0 && *g6mcok ){
    *trueSetFound = AreIdentical(g6mc, bestSet, ncll);      
    *trueSetFound = 0;
  }else{
    *trueSetFound = 0;
  }
  
  
  delete Rec;
  
} // END reconstruct3pi0_


/** my power function ***/
/*
double pow(double x, int n){
  double pow = 1;
  for(int i=0;i<n;i++){
    pow *= x;
  }
  return pow;
}
*/

/******************* distance between two points ********************/
double Dist(double * a, double * b){
  double d = 0;
  d += pow( a[0]-b[0] ,2.);
  d += pow( a[1]-b[1] ,2.);
  d += pow( a[2]-b[2] ,2.);
  d = sqrt( d );

  return d;
}

/**  A function to provide random downscaling capability  *******/
extern "C" {

  void downscale_(int * rate, int * result);

}

void downscale_(int * rate, int * result){
  //  srand(time(NULL));
  if(*rate==1){
    *result = 1;
    return;
  }
  double r = rand()/(double)RAND_MAX;
  if( r < 1./(*rate) ){
    *result = 1;
  }else{
    *result = 0;
  }
}




