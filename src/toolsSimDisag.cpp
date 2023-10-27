#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::NumericVector matrixSubcol(Rcpp::NumericMatrix x, int i0, int i1, int icol){ 
  
  // Determine the length of the output
  int n = i1-i0+1;
  
  // Create an output vector
  Rcpp::NumericVector out(n);
  
  // Loop through the matrix
  for(int z = 0; z < n; ++z) {
    out(z) = x(i0+z,icol);
  }
  
  return out;
}

// [[Rcpp::export]]
void rcpp_rprintf(NumericVector v){
  // printing values of all the elements of Rcpp vector  
  for(int i=0; i<v.length(); ++i){
    Rprintf("the value of v[%i] : %f \n", i, v[i]);
  }
}

// [[Rcpp::export]]
bool all_sug(LogicalVector x) {
  // Note the use of is_true to return a bool type.
  return is_true(all(x == TRUE));
}


// [[Rcpp::export]]
int find_row(NumericMatrix m, NumericVector v) {
  // a function that returns the index of the row which is equal to the vector vec
  for (int i = 0; i < m.nrow(); i++) {
    if (all_sug(m(i,_) == v)) {
      return i;
    }
  }
  return -1;
}

// [[Rcpp::export]]
NumericVector simPrecipOccurrences(int nLag, const NumericMatrix & matcomb, 
                                   NumericMatrix Qtrans, NumericVector rndNorm) {
  // a function that takes as input the number of Lags to consider, the length of simulated occurrences process, 
  // matcomb : a matrix with possible wet/dry combinations
  // Qtrans : a matrix n X nComb with quantile of each probability for each combination of rows in matcomb
  // rndNorm : a vector of random normal deviates that will be confronted to Qtranc_vec
  
  // define a vector of length nChainFit that will hold the occurrences process 
  int n = rndNorm.length();
  NumericVector Xt(n);
  
  // Fill the first time steps with random numbers (0 and 1 with proba 0.5)
  for (int t = 0; t < nLag; t++){
    Xt[t] = Rcpp::runif(1)[0]<0.5;
  }
  
  // confront each combination Qtranc_vec quatile to normal deviate
  NumericVector comb(nLag);
  
  for (int t = nLag; t < n; t++){
    // find nLag previous occurrences
    comb = Xt[Range(t-nLag,t-1)];
    
    // find which row of the matcomb matrix is equal to current combination
    int row = find_row(matcomb, comb);
    
    // confront to normal deviate, assign the value to the current occurrence 
    Xt[t] = rndNorm[t]<=Qtrans(t,row);
  }
  
  return Xt;
}


// [[Rcpp::export]]
NumericVector simPrecipOccurrences4Fitting(int nLag, int nChainFit, const NumericMatrix & matcomb, 
                                           NumericVector Qtrans_vec, NumericVector rndNorm) {
  // a function that takes as input the number of Lags to consider, the length of simulated occurrences process, 
  // matcomb : a matrix with possible wet/dry combinations
  // Qtrans_vec : a vector with quantile of each probability for each combination of rows in matcomb
  // rndNorm : a vector of random normal deviates that will be confronted to Qtranc_vec
  
  // define a vector of length nChainFit that will hold the occurrences process 
  int n = rndNorm.length();
  NumericVector Xt(n);
  
  // Fill the first time steps with random numbers (0 and 1 with proba 0.5)
  for (int t = 0; t < nLag; t++){
    Xt[t] = Rcpp::runif(1)[0]<0.5;
  }
  
  // confront each combination Qtranc_vec quatile to normal deviate
  NumericVector comb(nLag);
  
  for (int t = nLag; t < n; t++){
    // find nLag previous occurrences
    comb = Xt[Range(t-nLag,t-1)];
    
    // find which row of the matcomb matrix is equal to current combination
    int row = find_row(matcomb, comb);
    
    // confront to normal deviate, assign the value to the current occurrence 
    Xt[t] = rndNorm[t]<=Qtrans_vec[row];
  }
  
  return Xt[Range(n-nChainFit,n-1)];
}

// [[Rcpp::export]]
double pearsonrho(NumericVector x, NumericVector y) {
  // a function that return pearson rho
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
  double sum_x2 = 0.0, sum_y2 = 0.0;
  int n = x.length();
  
  for (int i = 0; i < n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
    sum_x2 += x[i] * x[i];
    sum_y2 += y[i] * y[i];
  }
  
  double numerator = n * sum_xy - sum_x * sum_y;
  double denominator = sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));
  
  return numerator / denominator;
}

// [[Rcpp::export]]
double corMarkovChain(NumericMatrix rndNorm,
                      NumericMatrix QtransMat, 
                      const NumericMatrix & matcomb, 
                      int nChainFit, int nLag){
  // initialise Xt with zeros
  int n = rndNorm.nrow();
  NumericMatrix Xt(n, 2);
  
  // sequence of previous nLag wet/dry states for one station and one time step
  NumericVector comb(nLag);
  
  // Fill the first time steps with random numbers (0 and 1 with proba 0.5)
  for (int t = 0; t < nLag; t++){
    Xt(t,0) = rndNorm(t,0)<0;
    Xt(t,1) = rndNorm(t,1)<0;
  }
  
  // Fill Xt
  for (int t = nLag; t < n; t++) {
    for (int st = 0; st < 2; st++) {
      // nLag past states
      comb = matrixSubcol(Xt,t-nLag,t-1,st);
      
      // find which row of the matcomb matrix is equal to current combination
      int iComb = find_row(matcomb, comb);
      
      // confront to normal deviate, assign the value to the current occurrence 
      Xt(t,st) = rndNorm(t,st)<=QtransMat(st,iComb);
    }
  }
  
  // correlation between number
  return pearsonrho(matrixSubcol(Xt,n-nChainFit,n-1,0),matrixSubcol(Xt,n-nChainFit,n-1,1));
}

// [[Rcpp::export]]
IntegerVector order_(NumericVector x) {
  NumericVector sorted = clone(x).sort();
  
  return match(sorted, x);
}

// [[Rcpp::export]]
bool any_sug(LogicalVector x){
  // Note the use of is_true to return a bool type
  return is_true(any(x == TRUE));
}

// [[Rcpp::export]]
NumericVector getrmsei(NumericMatrix Yobs, // Matrix of observations: nTobs*3 [days] x nStat [stations]
                       NumericMatrix Y3obs, // Matrix of observations amounts for 3-day periods: nTobs [3-days] x nStat [stations]
                       IntegerVector mObs, // index of the month (1..12) for the matrix of 3-day observations
                       IntegerVector cObs, // class (=1,2,3,4) of precip for observations
                       NumericMatrix Y3sim, // Matrix of simulated amounts for 3-day periods: nTsim [3-days] x nStat [stations]
                       IntegerVector mSim, // index of the month (1..12) for the matrix of 3-day simulations
                       IntegerVector cSim, // class (=1,2,3,4) of precip for simulations
                       int nLagScore, // lag=1
                       int i,
                       NumericMatrix Ysimilag){ // Ysim for the 3 days before the current i: 3xk){
  
  // Input/Output
  int nTobs = Y3obs.nrow();
  int nStat = Y3obs.ncol();
  
  // Locals
  const double& naVal=-9999.;
  double rmseIJ;
  NumericVector adimObs(nStat);
  NumericVector adimSim(nStat);
  NumericVector rmseI(nTobs);
  bool notSameClass;
  int jDay; // index of the observed day
  
  
  //======= Step 1: Minimize Score for this selection
  //sum of square differences for the closest fields on two time steps
  
  //for the first time steps
  if(i < nLagScore) {
    
    for (int j = 0; j < nTobs; j++) {
      
      //same month and class
      notSameClass = (mSim(i)/=mObs(j)) | (cSim(i)/=cObs(j));
      //if any observed value is missing in the observed prec. field or if the months do not correspond
      //we discard this observed day
      if((any_sug(Y3obs(j,_) == naVal)) | notSameClass) {
        rmseI(j) = 1E30;
      } else {
        //absolute differences between adimensioned precipitation for this day
        if(sum(Y3sim(i,_))==0) {
          adimSim = Y3sim(i,_);
        } else {
          adimSim = Y3sim(i,_)/sum(Y3sim(i,_));
        }
        
        if(sum(Y3obs(j,_))==0) {
          adimObs = Y3obs(j,_);
        } else {
          adimObs = Y3obs(j,_)/sum(Y3obs(j,_));
        }
        rmseI(j) = sum(abs(adimSim-adimObs));
      }
      
    }
  } else {
    
    //discard first elements
    for (int j = 0; j < nLagScore; j++) {
      rmseI(j) = 2E30;
    }
    
    
    //for the next days, compute score
    for (int j = nLagScore; j < nTobs; j++) {
      
      //same month and class
      notSameClass = (mSim(i)!=mObs(j)) | (cSim(i)!=cObs(j));
      //if any observed value is missing in the observed prec. field or if the months do not correspond
      //we discard this observed day
      if(any_sug(Y3obs(j,_) == naVal) | notSameClass) {
        rmseI(j) = 3E30;
      } else {
        //absolute differences between adimensioned precipitation for this day
        if(sum(Y3sim(i,_))==0) {
          adimSim = Y3sim(i,_);
        } else {
          adimSim = Y3sim(i,_)/sum(Y3sim(i,_));
        }
        
        if(sum(Y3obs(j,_))==0) {
          adimObs = Y3obs(j,_);
        } else {
          adimObs = Y3obs(j,_)/sum(Y3obs(j,_));
        }
        rmseIJ = sum(abs(adimSim-adimObs));
        
        
        // add differences for the previous days, just non na values
        // jDay is the index in the daily matrix just before the current
        // 3-day step. We  roll over the nLag previous daily steps starting
        // from jDay to jDay-nLag+1 (iLag = 0,...,(nLag-1))
        
        //index just before this 3-day period in the daily matrix
        jDay = j*3-1;
        
        for (int iLag = 0; iLag < nLagScore; iLag++) {
          if(any_sug(Yobs(jDay-iLag,_)==naVal)) {
            rmseIJ = 4E30;
            break;
          } else {
            if(sum(Ysimilag(2-iLag,_))==0) {
              adimSim = Ysimilag(2-iLag,_);
            } else {
              adimSim = Ysimilag(2-iLag,_)/sum(Ysimilag(2-iLag,_));
            }
            
            if(sum(Yobs(jDay-iLag,_))==0) {
              adimObs = Yobs(jDay-iLag,_);
            } else {
              adimObs = Yobs(jDay-iLag,_)/sum(Yobs(jDay-iLag,_));
            }
            
            rmseIJ = rmseIJ + sum(abs(adimSim-adimObs));
          }
        }
        
        rmseI(j) = rmseIJ;
      }
    }
  }
  return rmseI;
}


// [[Rcpp::export]]
List disag3DayGWexPrec(NumericMatrix Yobs, // Matrix of observations: nTobs*3 [days] x nStat [stations]
                       NumericMatrix Y3obs, // Matrix of observations amounts for 3-day periods: nTobs [3-days] x nStat [stations]
                       IntegerVector mObs, // index of the month (1..12) for the matrix of 3-day observations
                       IntegerVector cObs, // class (=1,2,3,4) of precip for observations
                       NumericMatrix Y3sim, // Matrix of simulated amounts for 3-day periods: nTsim [3-days] x nStat [stations]
                       IntegerVector mSim, // index of the month (1..12) for the matrix of 3-day simulations
                       IntegerVector cSim, // class (=1,2,3,4) of precip for simulations
                       int nLagScore){
  // a function that disagregate 3-day precipitations with daily precipitation
  // For each 3-day simulations at n gauges, we find the closest field in observations for which
  // 3-day precipitation structures are available.
  
  // Input/Output
  
  int nTobs = Y3obs.nrow();
  int nStat = Y3obs.ncol();
  int nTsim = Y3sim.nrow();
  
  // outputs
  NumericMatrix Ysim(nTsim*3,nStat); // Matrix of disag. simulated amounts: nTsim*3 [days] x nStat [stations]
  NumericMatrix codeDisag(nTsim,nStat); //Matrix indicating how it has been disaggregated
  
  // Locals
  const int& nBestField=10;
  const double& naVal=-9999.;
  double rmseIJ;
  NumericVector adimObs(nStat);
  NumericVector adimSim(nStat);
  NumericVector Yobs3D(3);
  NumericVector rmseI(nTobs);
  IntegerVector indBestRMSEI(nTobs);
  IntegerVector indBestFieldI(nBestField);
  bool notSameClass;
  int iDay; // index of the simulated day
  int jDay; // index of the observed day
  NumericMatrix Ysimilag(3,nStat);
  
  //For each simulated 3-day period
  for (int i = 0; i < nTsim; i++) {
    
    //index just before this 3-day period in the daily matrix
    iDay = i*3-1;
    
    // score
    if(i>0){
      for(int ilag = 0;ilag<3; ilag++){
        Ysimilag(ilag,_) = Ysim(iDay+ilag-2,_);
      }
      
    }
    rmseI = getrmsei(Yobs,Y3obs,mObs,cObs,Y3sim,mSim,cSim,nLagScore,i,Ysimilag);
    
    // get index of minimum RMSE values
    // !!! order_ starts from 1, which means that indBestFieldI is shifted by 1
    // in terms of index of the precipitation matrices
    indBestRMSEI = order_(rmseI)-1;
    indBestFieldI = indBestRMSEI[Range(0,nBestField-1)];
    
    //======= Step 3: Look at the different case and disaggregate =====
    for (int k = 0; k < nStat; k++) {
      
      //initialise code
      codeDisag(i,k) = naVal;
      
      
      //////////case 1: no occurrence for this period and this station, nothing to disaggregate.
      if(Y3sim(i,k)==0) {
        codeDisag(i,k) = -1000;
        for (int ii = iDay+1; ii <= (iDay+3); ii++) {
          Ysim(ii,k) = 0;
        }
        
        //////////case 2: loop at the closest fields, if there is the same 
        // number of occurrences, we take the observed
        //temporal structure (less restrictive than the exact same structure)
      } else {
        
        for (int j = 0; j < nBestField; j++) {
          
          //index just before this 3-day period in the daily matrix
          jDay = indBestFieldI(j)*3-1;
          
          //check if there is an observed precitation for the selected observed 
          // field and for this station
          double Y3obsjk = Y3obs(indBestFieldI(j),k);
          if(Y3obsjk>0) {
            
            //update code of disaggregation
            codeDisag(i,k) = j+1;
            //simulated precipitation for these 3 days are observed 
            // precipitation for the close field, rescaled
            //by the simulated cumulated value for these 3 days
            //for(int i3D = (jDay+1);i3D<=(jDay+3);i3D++){
            //  Yobs3D(i3D-jDay-1) = Yobs(i3D,k);
            //}
            Yobs3D = matrixSubcol(Yobs,jDay+1,jDay+3,k);
            
            // iDay is the index in the daily matrix just before the current
            // 3-day step
            for (int ii = iDay+1; ii <= (iDay+3); ii++) {
              Ysim(ii,k) = Yobs3D(ii-iDay-1)*Y3sim(i,k)/sum(Yobs3D);
            }
            
            //codes to analyse how large 3-day intensities are disaggregated
            if(any_sug(matrixSubcol(Ysim,iDay+1,iDay+3,k)>max(Yobs(_,k)))){
              codeDisag(i,k) = codeDisag(i,k) + 10000;
            }
            //get out of the loop "for (int j = 0; j < nBestField; j++)"
            break;
          }
        }
      }
      
      //////////case 3: if we did not find similar structure then we find, for
      // this station, days with similar amounts
      if(codeDisag(i,k)==naVal) {
        //for the first time step
        if(i < nLagScore) {
          for (int j = 0; j < nTobs; j++) {
            //same month and class
            notSameClass = (mSim(i)/=mObs(j)) | (cSim(i)/=cObs(j));
            //if any observed value is missing in the observed prec. field or if the months do not correspond
            //we discard this observed day
            if((Y3obs(j,k) == naVal) | notSameClass) {
              rmseI(j) = 1E30;
            } else {
              rmseI(j) = abs(Y3sim(i,k)-Y3obs(j,k));
            }
          }
        } else {
          //discard first elements
          for (int j = 0; j < nLagScore; j++) {
            rmseI(j) = 1E30;
          }
          
          //for the next days, compute score
          for (int j = nLagScore; j < nTobs; j++) {
            //same month and class
            notSameClass = (mSim(i)!=mObs(j)) | (cSim(i)!=cObs(j));
            //if any observed value is missing in the observed prec. field or if the months do not correspond
            //or if there is no observed precip for this day
            //we discard this observed day
            if(((Y3obs(j,k) == naVal) | (Y3obs(j,k) == 0)) | notSameClass) {
              rmseI(j) = 1E30;
            } else {
              rmseIJ = abs(Y3sim(i,k)-Y3obs(j,k));
              //add differences for the previous days, just non na values
              //index just before this 3-day period in the daily matrix
              jDay = j*3-1;
              if(any_sug(matrixSubcol(Yobs,jDay-nLagScore+1,jDay,k)==naVal)) {
                rmseIJ = 1E30;
              } else {
                for (int iLag = 0; iLag < nLagScore; iLag++) {
                  rmseIJ = rmseIJ + abs(Ysim(iDay-iLag,k)-Yobs(jDay-iLag,k));
                }
              }
              rmseI(j) = rmseIJ;
            }
          }
        }
        
        // get index of minimum RMSE values
        // !!! order_ starts from 1, which means that indBestFieldI is shifted by 1
        // in terms of index of the precipitation matrices
        indBestRMSEI = order_(rmseI)-1;
        indBestFieldI = indBestRMSEI[Range(0,nBestField-1)];
        
        codeDisag(i,k) = 2000;
        int rndDay = floor(Rcpp::runif(1)[0]*10); 
        int j3Day = indBestFieldI(rndDay);
        jDay = j3Day*3-1;
        Yobs3D = matrixSubcol(Yobs,jDay+1,jDay+3,k);
        for(int ii = iDay+1;ii<=(iDay+3);ii++){
          Ysim(ii,k) = Yobs3D(ii-iDay-1)*Y3sim(i,k)/sum(Yobs3D);
        }
        
        //codes to analyse how large 3-day intensities are disaggregated
        if(any_sug(matrixSubcol(Ysim,iDay+1,iDay+3,k)>max(Yobs(_,k)))){
          codeDisag(i,k) = codeDisag(i,k) + 10000;
        }
      }
    }
    
    
  }
  
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("codeDisag") = codeDisag,
                                      Rcpp::Named("Ysim") = Ysim);
  return res;
}
