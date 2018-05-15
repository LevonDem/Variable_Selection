
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ C++ Code for overlapping group variable selection ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
double pi = 3.141592653589793238462643383280;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~~ Sample from the univariate normal distribution  ~~
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
double rnormalCS(Function f, double x1, double x2) {
  double out = 0.0;
  List temp;
  
  temp = f(1, x1, x2);
  out = as<double>(temp[0]);
  
  return out;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Sample from the inverse gamma distribution  ~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
double rigammaCS(Function f, double x1, double x2) {
  double out = 0.0;
  List temp;
  
  temp = f(1, x1, x2);
  out = as<double>(temp[0]);
  
  return out;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Sample from the beta distribution  ~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
double rbetaCS(Function f, double x1, double x2) {
  double out = 0.0;
  List temp;
  
  temp = f(1, x1, x2);
  out = as<double>(temp[0]);
  
  return out;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Compute pdf of multivariate normal distribution, log scale ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
double logmvnpdfCS(const mat x, const mat mu, const mat Sigma, const int k) {
  double t1 = 0.0;
  double t2 = 0.0;
  mat invSigma = inv(Sigma);
  mat A = -0.5 * trans(x - mu) * invSigma * (x - mu);
  
  //t1 = -0.5 * k * log(2.0 * pi) - 0.5 * log(det(Sigma));
  t1 = - 0.5 * log(det(Sigma));
  t2 = A[0];
  
  return t1 + t2;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Binary sampling: Prob(X=1) ~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
int sampleBinaryS(double prob){
                
  double coinFlip = R::runif(0,1);
  int output = 0;
  
  if(coinFlip <= prob){
    output = 1;
  }
      
  return output;    
  
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Sample from a vector ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
int sampleCS(IntegerVector x){
                     
  NumericVector prob = NumericVector::create();                   
  IntegerVector ret = Rcpp::RcppArmadillo::sample(x, 1, false, prob);
  
  return as<int>(ret);
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Set difference (vector version) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
vec setdiffCS(IntegerVector x, int y){   
  int j = 0;
  int n = x.size();
  vec out(n-1);
  
  for(int i = 0; i < n; ++i){
    if(x[i] != y){  
      out[j] = x[i];
      j++;
    }
  }
  return out;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Set difference (matrix version) ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
mat setdiffC2S(IntegerVector x){ 
  int n = x.size();
  mat y(n, n - 1);

  for(int i = 0; i < n; ++i){
    y.row(i) = setdiffCS(x, x[i]).t(); 
  }
  
  return y;
  
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~~~ Find group overlap structure ~~~
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
List overlap(List varInd){
  
  /* Find the overlap structure and size of the overlap 
  structure of the groups and variables */
  int nGroups = varInd.size();
  List overSize(nGroups), overStruct(nGroups);
  vec diffvec(nGroups);
  CharacterVector x(nGroups), y(nGroups), z(nGroups);
  IntegerVector w(nGroups), v(nGroups);
  int i = 0, j = 0, k = 0;
  for(i = 0; i < nGroups; i++){
    overStruct[i] = 0;
    overSize[i]   = 0;
    diffvec    = setdiffCS(seq_len(nGroups), (i + 1));
    for(j = 0; j < (nGroups - 1); j ++){
      k = diffvec[j] - 1;
      x = as<CharacterVector>(varInd[i]);
      y = as<CharacterVector>(varInd[k]);
      z = intersect(x, y);
      if(z.size() > 0){
        w = as<IntegerVector>(overStruct[i]);
        v = as<IntegerVector>(overSize[i]);
        w.push_back(k+1);
        v.push_back(z.size());
        overStruct[i] = w;
        overSize[i] = v;
      }
      
      w = as<IntegerVector>(overStruct[i]);
      v = as<IntegerVector>(overSize[i]);
      
      overStruct[i] = w[w != 0];
      overSize[i]   = v[v != 0];
      
      w = as<IntegerVector>(overStruct[i]);
      v = as<IntegerVector>(overSize[i]);
      
      if(w.size() == 0) overStruct[i] = 0;
      if(v.size() == 0) overSize[i] = 0;
      
    }
  }
  
  return List::create(
    Named("struct") = overStruct, 
    Named("size") = overSize
  );
  
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Adding elements to a list ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* This function takes in a list and outputs a list, with gamma added to list element n 
NOTE: this function directly alters the input x, DONT ASSIGN IT TO A VARIABLE (just call it)*/
// [[Rcpp::export]] 
List resizeListS(List x, int n, int gamma){

    NumericVector z = x[n];
    NumericVector q(z.size() + 1);
    
    for(int i = 0; i < (q.size() - 1); i++){
      q[i] = z[i];
    }

    q[q.size() - 1] = gamma;

    x[n] = q;
    return x;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Adding elements to a list ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* This function takes in a list and outputs a list, with gamma added to list element n 
NOTE: this function directly alters the input x, DONT ASSIGN IT TO A VARIABLE (just call it)*/
// [[Rcpp::export]] 
List resizeListDoubleS(List x, int n, double beta){

    NumericVector z = x[n];
    NumericVector q(z.size() + 1);
    
    for(int i = 0; i < (q.size() - 1); i++){
      q[i] = z[i];
    }

    q[q.size() - 1] = beta;

    x[n] = q;
    return x;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Extract last entry from each list element ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]] 
IntegerVector getRecentS(const List x){

    int n = x.size();
    IntegerVector y(n);
    IntegerVector z = x[0];

    for(int i = 0; i < n; i++){
      z = x[i];
      y[i] = z[z.size() - 1];
    }
    
    return y;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Extract last entry from each list element ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]] 
NumericVector getRecentDoubleS(const List x){

    int n = x.size();
    NumericVector y(n);
    NumericVector z = x[0];
    
    for(int i = 0; i < n; i++){
      z = x[i];
      y[i] = z[z.size() - 1];
    }
  
    return y;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Extract last entry from a Numeric Vector ~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]] 
double getRecentVectorS(const NumericVector x){

    int    n = x.size();
    double y = 1.0;

    if(n == 1){
      y = x[0];
    }
    else{
      y = x[n-1];
    }
    
    return y;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Markov Random Field prior on group indicators ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]] 
double groupPriorS(const IntegerVector etaTilde, 
                   const IntegerVector overlappingGroups, 
                  const IntegerVector sizeoverlapGroups, 
                  const int MRF, 
                  const double MRFweight){
                           
  int j = 0;
  int currentEta = 0;
  int otherEtaSum = 0;  
  double thetaL = 0.50;
  double groupEffect = 0.0;
  
  if(MRF == 1){
    if(overlappingGroups[0] > 0){
      for(j = 0; j < overlappingGroups.size(); j++){
        //#vec currentEtaVec = etaTilde[overlappingGroups[j] - 1];
        //#currentEta = currentEtaVec[currentEtaVec.size() - 1];
        currentEta = etaTilde[overlappingGroups[j] - 1];
        otherEtaSum += currentEta;
        groupEffect += MRFweight * sizeoverlapGroups[j] * (2 * currentEta - 1);
        //groupEffect  += (2 * currentEta - 1);
      }
    }
    thetaL = exp(-groupEffect) / (exp(-groupEffect) + exp(groupEffect)); // P(eta = 0)
  }
  
  return thetaL; 
    
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Sample sigma from the half-Cauchy distribution  ~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
double sampleSigmaCauchyS(Function rigamma, const int n, double q, const double S, const double oldSigma) {
  
  int tick = 0;
  double coinFlip = 0.0;
  double out = oldSigma;
  double sigmaProposed = 1.0;

  if(q <= pow(10, -20)){
	q = pow(10, -20);
  }
  
  // Implement rejection sampling
  while(tick < 100){
    tick += 1; 
    sigmaProposed = rigammaCS(rigamma, (n - 1.0)/2.0, q / 2.0); // this is sigma SQUARED
    //sigmaProposed = rigammaCS(rigamma, (n - 2.0)/2.0, q / 2.0); // this is sigma SQUARED
    coinFlip = R::runif(0,1);
    if(coinFlip <= (pow(S,2) / (pow(S,2) + sigmaProposed))){
    //if(coinFlip <= (pow(S,2) / (pow(S,2) + pow(sigmaProposed, 2)))){
      out = sqrt(sigmaProposed);
      break;
    }
  }
  
  return out;
  
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Sample sigma0 and sigma1 from the half-Cauchy distribution  ~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]]
double sampleSigma01CauchyS(Function rigamma, const int N, const double S, double betaSq, const double oldSigma) {
  
  int tick = 0;
  double coinFlip = 0.0;
  double sigmaProposed = 1.0;
  double newSigma = oldSigma;

  if(betaSq <= pow(10, -20)){
	  betaSq = pow(10, -20);
  }
  
  if(N == 0){
    coinFlip = R::runif(0,1);
    newSigma = S * tan(pi * coinFlip / 2.0);
  }
  else{
    while(tick < 100){
      tick += 1; 
      sigmaProposed = rigammaCS(rigamma, 0.5 * (N + 1), 0.5 * betaSq);
      coinFlip = R::runif(0,1);
      if(coinFlip <= pow(S, 2) / (sigmaProposed + pow(S, 2))){
        newSigma = sqrt(sigmaProposed);
        break;
      }
    }
  }
  
  return newSigma;
  
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Use sampling to approximate the log of the normalizing constant ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export]] 
double normalizingConstantUniformS(const int numSim, const int pl, const double thetaL, 
                           const double sigmaTilde, const vec s0, const vec s1, 
                           const IntegerVector gammaCW, const mat XGroup, const mat Rl, 
                           const mat Bl){
  /* Group size is fixed here, so we can declare everything outside our loop*/
  int sim = 0;
  int i = 0;
  int j = 0;
  int l = 0;
  double K = 1.0;
  double gammaSum = 0.0;
  double sampleMean = 0.0;
  double t = 0.0;  
  double logZfl = 1.0;
  mat b;
  mat Al;
  mat Sl;
  mat Vl;
  IntegerVector sampledGammas = gammaCW;
  vec zeroVec;
    
  if(pl > 10){  
    zeroVec.zeros(numSim);
    for(sim = 0; sim < numSim; sim++){
      
      for(i = 0; i < pl; i++){
        sampledGammas[i] = sampleBinaryS(0.5);
      }
          
      Vl = eye(pl, pl); 
      gammaSum = 0.0;
      for(l = 0; l < pl; l++){
        //sampledGammas[l] = sampleBinaryS(0.5);
        Vl(l,l) = 1.0 / (pow(s0[l], 2 * (1 - sampledGammas[l])) * pow(s1[l], 2 * sampledGammas[l]));
        gammaSum += sampledGammas[l];
      }
      Al = Vl + trans(XGroup) * XGroup / pow(sigmaTilde, 2);  
      Sl = trans(Rl) * XGroup * inv(Al) * trans(XGroup) * Rl;      
      t = pow(thetaL, gammaSum) * pow(1.0 - thetaL, pl - gammaSum);
      
      zeroVec[sim] = 0.5 * (log(det(Vl)) - log(det(Al))) + Sl[0]/(2*pow(sigmaTilde,4)) + log(t); 
    }
    
  }
  else{
    zeroVec.zeros(pow(2, pl));
    for(j = 0; j < pow(2, pl); j++){ 
      b = Bl.row(j);
      Vl = eye(pl, pl); 
      gammaSum = 0;
      for(l = 0; l < pl; l++){
        Vl(l,l) = pow(s0[l], 2 * (1 - b[l])) * pow(s1[l], 2 * b[l]);
        gammaSum += b[l];
      }      
      Al = inv(inv(Vl) + trans(XGroup) * XGroup / pow(sigmaTilde, 2));
      Sl = trans(Rl) * XGroup * Al * trans(XGroup) * Rl; 
      t = pow(thetaL, gammaSum) * pow(1.0 - thetaL, pl - gammaSum);
    
      zeroVec[j] = 0.5 * (log(det(Al)) - log(det(Vl))) + Sl[0]/(2*pow(sigmaTilde,4)) + log(t);     
    
    }    
  }
  
  K = max(zeroVec);
  if(pl > 10){     
    // Subtract each term by max    
    for(sim = 0; sim < numSim; sim++){
      zeroVec[sim] = zeroVec[sim] - K;
      sampleMean += exp(zeroVec[sim]);
    }
    logZfl = pl * log(2) - log(numSim) + K + log(sampleMean);  
  }
  if(pl <= 10){     
    // Subtract each term by max    
    for(j = 0; j < pow(2, pl); j++){
      zeroVec[j] = zeroVec[j] - K;
      sampleMean += exp(zeroVec[j]);
    }
    logZfl = K + log(sampleMean);  
  }
   
  return logZfl;
  
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ Component-wise Gibbs Sampler ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* This function agrees with a similar version I programmed in R*/
// [[Rcpp::export]]              
List CWGibbsS(const mat Xgroup, const mat Rl, const vec betaGroup, const IntegerVector gammaGroup,
               const int pl, const vec sigma0Group, const vec sigma1Group, const double sigmaTilde, 
               const int CWiter, const int r, const int s,
               Function rigamma, Function rnorm, const double S) {
  
  /* Group size is fixed here, so we can declare everything outside our loop*/                        
  double tauTilde = 0.0;
  double gammaPrior = 0.0;
  double sigma0 = 1.0; 
  double sigma1 = 1.0; 
  double denominator = 0.0;
  double sigma = sigmaTilde;
  double sigmaStar = 0.0;
  double rli = 0.0;  
  double rl0 = 0.0;
  double rl1 = 0.0;
  double den0 = 0.0;
  double den1 = 0.0;
  double logzli = 0.0;
  double zli = 0.0;
  double pGamma = 0.0;
  int gammaHat = 0;  
  int otherGammaSum = 0; 
  int n = Xgroup.n_rows;
  IntegerVector currentGammaVec = gammaGroup;
  IntegerVector M1 = seq_len(pl);
  mat XX;
  mat RX;
  mat Rli = Rl;  
  mat Res;
  mat resInnerProd;
  vec betaNew = betaGroup;
  vec otherGammas;

  // Begin CGS
  for(int k = 0; k < CWiter; k++){
    
    for(int j = 0; j < pl; j++){
      
      // Save the sum of the other gammas in this group
      otherGammaSum = 0;
      for(int i = 0; i < pl; i++){
        otherGammaSum += currentGammaVec[i];
      }
      otherGammaSum = otherGammaSum - currentGammaVec[j];
      
      // Treat Rl as new Y in this group and compute residual
      Rli = Rl;
      for(int i = 0; i < pl; i++){
        Rli = Rli - (Xgroup.col(i) * betaNew[i]); 
      }
      Rli = Rli + Xgroup.col(j) * betaNew[j]; 

      // Extract the sd's for each group
      sigma0 = sigma0Group[j];
      sigma1 = sigma1Group[j];
    
      // Compute X'X and R'X 
      XX = trans(Xgroup.col(j)) * Xgroup.col(j);
      RX = trans(Rli) * Xgroup.col(j);

      // Compute r0 and r1
      den0 = pow(sigma0, 2) * XX[0] + pow(sigma, 2);
      den1 = pow(sigma1, 2) * XX[0] + pow(sigma, 2);
    
      rl0 = pow(sigma0, 2) * pow(RX[0], 2) / (den0);
      rl1 = pow(sigma1, 2) * pow(RX[0], 2) / (den1);
    
      // Finally, compute zli
      logzli = 0.5 * (log(den0) - log(den1)) + (rl1 - rl0) / (2 * pow(sigma, 2));
      zli    = exp(logzli);

      // For numerical insurance:
      if(zli >= pow(10, 200))
        zli = pow(10, 200);
      
      // Sample gamma conditional on eta=1 and beta
      gammaPrior = (s + pl - otherGammaSum - 1.0) / (r + s + pl - 1.0); // P of gamma = 0
      pGamma     = ((1 - gammaPrior) * zli) / (((1 - gammaPrior) * zli) + gammaPrior);
      gammaHat   = sampleBinaryS(pGamma);

      // Sample beta conditional on eta=1 and gamma
      tauTilde    = (1 - gammaHat) * sigma0 + gammaHat * sigma1; 
      denominator = XX[0] * pow(tauTilde, 2) + pow(sigma, 2); 
      sigmaStar   = sigma * tauTilde / sqrt(denominator); 
      rli         = (RX[0] * pow(tauTilde, 2)) / denominator; 
    
      // Save new values of gamma and beta
      betaNew[j]         = R::rnorm(rli, sigmaStar);
      currentGammaVec[j] = gammaHat; 

    }

    // Compute residual vector
    Res = Rl - (Xgroup * betaNew);
    resInnerProd = (trans(Res) * Res);
    
    // Sample sigma from posterior distribution
    sigma = sampleSigmaCauchyS(rigamma, n, resInnerProd(0), S, sigma);
  
  }
  
  //Rcout << "resInnerProd" << resInnerProd(0) << std::endl;
  
  // Return 1 sample of gamma vector and beta vector for this group
  return List::create(
    //Named("CWgammas") = gammaNew, 
    Named("CWgammas") = currentGammaVec,
    Named("CWbetas")  = betaNew
  );


}

  
// Rcout << "betaNew: " << betaNew << std::endl;
// Rcout << "gammaHat: " << gammaHat << std::endl;


// [[Rcpp::export]]                 
List SGW(const List Xtilde, const mat Y, List betaTilde, IntegerVector etaTilde, List gammaTilde, 
                 const IntegerVector p, List sigma0, List sigma1, 
                 const int CWiter, const double sigmaTilde, const List overlapStructure, 
                 const List sizeOverlapStructure, const int numSim, const double BurnIn, 
                 NumericVector theta, const int MRF, const double MRFweight, const List B,
                 Function rbeta, Function rigamma, Function rnorm, const double S0, const double S1, const double S,
                 const vec r, const vec s, IntegerVector etaCount) { // BurnIn should be a percentage
                         
  int nGroups = p.size();
  int currentGamma = 0;
  int gammaSum = 0;
  int etaHat = 0;
  int curEtaInt = 0;
  int gammaTotal = 0;
  int i = 0;
  int j = 0;
  int l = 0;
  int N0 = 0;
  int N1 = 0;
  double Alhat = 0.0;
  double Dlhat = 0.0;
  double betaSq0 = 0.0;
  double betaSq1 = 0.0;
  double currentBeta = 0.0;
  double newSigma0 = 1.0;
  double newSigma1 = 1.0;
  double piece1 = 0.0;
  double piece2 = 0.0;
  double piece3 = 0.0;
  double thetaL = 0.0;
  double T1 = 0.0;
  double T2 = 0.0;
  double logZfl = 1.0;
  IntegerVector XtildeLengthSeq = seq_len(p.size());
  IntegerVector currentEtaVec;
  IntegerVector curEta;
  
  for(i = 0; i < nGroups; i++){
    
    // Save the indices of the other groups
    vec otherGroups = setdiffCS(XtildeLengthSeq, i + 1) - 1;
    mat XGroup = Xtilde[i];
    mat Bl = B[i];
    NumericVector betaGroup = getRecentDoubleS(betaTilde[i]);
    IntegerVector gammaGroup = gammaTilde[i];
    
    // Compute residual vector
    mat Rl = Y; 
    for(j = 0; j < (nGroups - 1); j++){
      Rl = Rl - as<mat>(Xtilde[otherGroups[j]]) * as<vec>(getRecentDoubleS(betaTilde[otherGroups[j]]));
    }
    
    // Save group indicator and group variances
    curEtaInt = etaTilde[i];
    vec s0 = as<vec>(sigma0[i]);
    vec s1 = as<vec>(sigma1[i]);
    
    // Save gamma and beta for this group
    IntegerVector gammaCW = gammaGroup;
    vec betaCW = betaGroup;
    
    if(curEtaInt == 0){
      
      // If current indicator of group is 0, sample beta and gamma via the 
      // CWGS, then decide whether or not to accept this proposal
      
      List CW = CWGibbsS(XGroup, Rl, betaCW, gammaGroup, p[i], s0, s1, sigmaTilde, CWiter, r[i], s[i], rigamma, rnorm, S);
      IntegerVector gammaCW = (CW[0]);
      betaCW = as<vec>(CW[1]);
      
      // If all proposed gammas are 0, zero out the group
      if(sum(gammaCW) == 0){
        etaTilde[i] = 0;
        gammaGroup  = 0 * gammaGroup;
        for(int j = 0; j < p[i]; j++)
          resizeListDoubleS(betaTilde[i], j, 0);
      }
      else{
        
        // If at least 1 gamma is not 0, propose to set eta = 1
        logZfl = normalizingConstantUniformS(numSim, p[i], theta[i], sigmaTilde, s0, s1, gammaCW, XGroup, Rl, Bl);
        mat Vl = eye(p[i], p[i]); 
        for(l = 0; l < p[i]; l++){
          Vl(l,l) = pow(s0[l], 2 * (1 - gammaCW[l])) * pow(s1[l], 2 * gammaCW[l]);
          gammaSum += gammaCW[l];
        }
        mat Al = inv(Vl) + trans(XGroup) * XGroup / pow(sigmaTilde, 2); 
        mat invAl = inv(Al);
        mat Sl = trans(Rl) * XGroup * invAl * trans(XGroup) * Rl;    
        mat mu = invAl * trans(XGroup) * Rl / pow(sigmaTilde, 2);
        vec zeroVec;
        zeroVec.zeros(p[i]); 
        
        //Rcout << "logZfl " << logZfl << std::endl; //

        piece1 = -0.5*(log(det(Al)) + log(det(Vl))) + (Sl[0]/(2*pow(sigmaTilde,4))) - logZfl;
        piece2 = logmvnpdfCS(betaCW, mu, invAl, p[i]);
        piece3 = logmvnpdfCS(betaCW, zeroVec, Vl, p[i]); 
        
        mat Cl = trans(betaCW) * trans(XGroup) * XGroup * betaCW - trans(Rl) * XGroup * betaCW - trans(betaCW) * trans(XGroup) * Rl;
        
        thetaL = groupPriorS(etaTilde, as<IntegerVector>(overlapStructure[i]), as<IntegerVector>(sizeOverlapStructure[i]), MRF, MRFweight);     
        //thetaL = 0.5;
        
        T1 = log(((1 - thetaL) / (thetaL))) + (-1/(2*pow(sigmaTilde, 2))) * Cl[0];
        T2 = piece3 - piece1 - piece2;
        
        // Compute log probability of group activation
        Alhat = T1 + T2;
        //Alhat = 1.0; // REMOVE
        
        // Decide whether or not to keep proposal
        if(Alhat >= 0){
          etaHat = 1;
        }
        else{
          etaHat = sampleBinaryS(exp(Alhat));
        }
        
        // Save variables
        if(etaHat == 1){
          etaTilde[i] = 1;
          gammaGroup  = gammaCW;
          for(int j = 0; j < p[i]; j++)
            resizeListDoubleS(betaTilde[i], j, betaCW[j]);
        }
        else{
          etaTilde[i] = 0;
          gammaGroup  = 0 * gammaGroup;
          for(int j = 0; j < p[i]; j++)
            resizeListDoubleS(betaTilde[i], j, 0);  
        }
      }
    } // End of if(curEtaInt == 0)
    else{
      
      // If current indicator of group is 1, let beta and gamma = 0
      // and decide whether or not to accept this proposal
      
      // If all current gammas are 0, zero out the group
      if(sum(gammaCW) == 0){
        etaTilde[i] = 0;
        gammaGroup  = 0 * gammaGroup;
        for(int j = 0; j < p[i]; j++)
          resizeListDoubleS(betaTilde[i], j, 0);
      }
      else{
      
        logZfl = normalizingConstantUniformS(numSim, p[i], theta[i], sigmaTilde, s0, s1, gammaCW, XGroup, Rl, Bl);
        
        mat Vl = eye(p[i], p[i]); 
        for(l = 0; l < p[i]; l++){
          Vl(l,l) = pow(s0[l], 2 * (1 - gammaCW[l])) * pow(s1[l], 2 * gammaCW[l]);
          gammaSum += gammaCW[l];
        }
        mat Al = inv(Vl) + trans(XGroup) * XGroup / pow(sigmaTilde, 2); 
        mat invAl = inv(Al);
        mat Sl = trans(Rl) * XGroup * invAl * trans(XGroup) * Rl;    
        mat mu = invAl * trans(XGroup) * Rl / pow(sigmaTilde, 2);
        vec zeroVec;
        zeroVec.zeros(p[i]); 

        piece1 = -0.5*(log(det(Al)) + log(det(Vl))) + (Sl[0]/(2*pow(sigmaTilde,4))) - logZfl;
        piece2 = logmvnpdfCS(betaCW, mu, invAl, p[i]);
        piece3 = logmvnpdfCS(betaCW, zeroVec, Vl, p[i]); 
        
        mat Cl = trans(betaCW) * trans(XGroup) * XGroup * betaCW - trans(Rl) * XGroup * betaCW - trans(betaCW) * trans(XGroup) * Rl;
        
        thetaL = groupPriorS(etaTilde, as<IntegerVector>(overlapStructure[i]), as<IntegerVector>(sizeOverlapStructure[i]), MRF, MRFweight);     
        //thetaL = 0.5;
        
        T1 = log(((1 - thetaL) / (thetaL))) + (-1/(2*pow(sigmaTilde, 2))) * Cl[0];
        T2 = piece3 - piece1 - piece2;
        
        // Compute log probability of group activation
        Alhat = T1 + T2;
        Dlhat = -Alhat;
        
        // Decide whether or not to keep proposal
        if(Dlhat >= 0){
          etaHat = 0;
        }
        else{
          etaHat = 1 - sampleBinaryS(exp(Dlhat));
        }
        
        // Save variables
        if(etaHat == 0){
          etaTilde[i] = 0;
          gammaGroup  = 0 * gammaCW;
          for(int j = 0; j < p[i]; j++)
            resizeListDoubleS(betaTilde[i], j, 0);
        }
        else{
          
          // Resample gamma and beta using CWGS
          List CW = CWGibbsS(XGroup, Rl, betaCW, gammaGroup, p[i], s0, s1, sigmaTilde, CWiter, r[i], s[i], rigamma, rnorm, S);
          IntegerVector gammaCW = (CW[0]);
          betaCW = as<vec>(CW[1]);
          
          // If all proposed gammas are 0, zero out the group
          if(sum(gammaCW) == 0){
            etaTilde[i] = 0;
            gammaGroup  = 0 * gammaGroup;
            for(int j = 0; j < p[i]; j++)
              resizeListDoubleS(betaTilde[i], j, 0);
          }
          else{
            etaTilde[i] = 1;
            gammaGroup  = gammaCW;
            for(int j = 0; j < p[i]; j++)
              resizeListDoubleS(betaTilde[i], j, betaCW[j]);
          }
        }
      }
    }
    
    // Save eta count
    etaCount[i] = etaTilde[i];
    //Rcout << "gammaCW " << gammaCW << std::endl;
    
  }
  
  return List::create(
    Named("etaNew") = etaTilde, 
    Named("gammaNew") = gammaTilde
    //Named("betaNew") = betaTilde
  );


}



// Rcout << "betaGroup: " << betaGroup[0] << std::endl;


// [[Rcpp::export]]                 
List outerLoopCS(const List Xtilde, const mat Y, List betaTilde, 
                List variableIndicator, List params, 
                Function rbeta, Function rigamma, Function rnorm,
                vec r, vec s){ 

  // Define model parameters
  int CWiter       = params[0];
  int outerIter    = params[1];
  int numSim       = params[2];
  int MRF          = params[3];
  double BurnIn    = as<double>(params[4]);
  double S         = as<double>(params[5]);
  double S0        = as<double>(params[6]);
  double S1        = as<double>(params[7]);
  double MRFweight = as<double>(params[8]);
  List B           = params[9];
  
  int i = 0;
  int j = 0;
  int k = 0;
  int n = Y.n_rows;                
  int nGroups = Xtilde.size();
  int currentEta = 0;
  int currentGamma = 0;
  int gammaTotal = 0;
  int N0 = 1;
  int N1 = 1;
  int p0 = 1; //
  double betaSq0 = 0.0;
  double betaSq1 = 0.0;
  double currentBeta = 1.0;
  double currentSigma = 1.0; 
  double newSigma = 1.0;
  double newSigma0 = 1.0;
  double newSigma1 = 1.0;
  //double p0 = 1.0;
  double validLoops = 1.0;
  IntegerVector g1;
  IntegerVector g2;
  mat Res;
  mat resInnerProd;
  NumericVector etaFinal(nGroups, 0.0);
  NumericVector sigmaTildeVec(1, 1.0);
  //NumericVector sigmaTildeVec = 1.0;
  NumericVector theta(nGroups, 0.5);
  IntegerVector etaCount(nGroups, 0);
  IntegerVector etaCumulative(nGroups, 0);
  IntegerVector etaTilde(nGroups, 1);
  IntegerVector p(nGroups);
  List gammaTilde(nGroups);
  List gammaCumulative(nGroups);
  List gammaFinal(nGroups);
  List groupInfo = overlap(variableIndicator);
  List overlapStructure = groupInfo[0];
  List sizeOverlapStructure = groupInfo[1];
  List sigma0(nGroups);
  List sigma1(nGroups);

  // Initialize lists for group and variable indicators
  for(i = 0; i < nGroups; i++){
    CharacterVector curVar = variableIndicator[i];
    p[i] = curVar.size();
    p0 = p[i];
    IntegerVector z(p0, 1);
    etaTilde[i] = 1;
    gammaTilde[i] = z;
    gammaCumulative[i] = z;
  }
  
  // Initialize lists for variances
  for(j = 0; j < nGroups; j++){
    mat varMat(p[j], 1);
    varMat.fill(1.0);
    sigma0[j] = 0.01 * varMat;
    sigma1[j] = 3.00 * varMat;
  }
  
  for(i = 0; i < outerIter; i++){
    Rcout << "Outer loop #" << (i + 1) << std::endl;
    currentSigma = getRecentVectorS(sigmaTildeVec);

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Perform sample-based variable selection */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    SGW(Xtilde, Y, betaTilde, etaTilde, gammaTilde, p, sigma0, sigma1, 
                CWiter, currentSigma, overlapStructure, sizeOverlapStructure, numSim, 
                BurnIn, theta, MRF, MRFweight, B, rbeta, rigamma, rnorm, S0, S1, S, r, s, etaCount);

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Sample theta from posterior distribution */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    
    for(k = 0; k < nGroups; k++){
      gammaTotal = 0;
      currentEta = etaTilde[k];
      IntegerVector gammaGroup = gammaTilde[k];
      for(j = 0; j < p[k]; j++){
        currentGamma = gammaGroup[j];
        gammaTotal += currentGamma;
      }
      
      theta[k]   = rbetaCS(rbeta, currentEta * gammaTotal + r[k], currentEta * (p[k] - gammaTotal) + s[k]);
    }
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Sample component variances from posterior distributions */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    
    for(k = 0; k < nGroups; k++){
      betaSq0    = 0.0;
      betaSq1    = 0.0;
      gammaTotal = 0;
      currentEta = etaTilde[k];
      IntegerVector gammaGroup = gammaTilde[k];
      NumericVector finalBeta  = getRecentDoubleS(betaTilde[k]);
      for(j = 0; j < p[k]; j++){
        currentBeta  = finalBeta[j];
        currentGamma = gammaGroup[j];
        gammaTotal += currentGamma;
        betaSq0 += pow(currentBeta, 2) * (1.0 - currentGamma);
        betaSq1 += pow(currentBeta, 2) * (currentGamma);
      }
      N1 = gammaTotal;
      N0 = p[k] - N1;
      
      mat curSigma0 = sigma0[k];
      mat curSigma1 = sigma1[k];
    
      // Sample from Half-Cauchy
      newSigma0 = sampleSigma01CauchyS(rigamma, currentEta * N0, S0, betaSq0, curSigma0[0]);
      newSigma1 = sampleSigma01CauchyS(rigamma, currentEta * N1, S1, betaSq1, curSigma1[0]);
      //newSigma0 = sqrt(rigammaCS(rigamma,  (N0 + 10.0) / 2.0, (betaSq0 + 1.0) /2.0));
      //newSigma1 = sqrt(rigammaCS(rigamma,  (N1 + 1.0) / 2.0, (betaSq1 + 10.0) /2.0));
  
      // Constrain sigma0 to be less than sigma1
      if(newSigma0 > 0.1 * newSigma1)
        newSigma0 = 0.1 * newSigma1;

      // Save values
      mat varMat(p[k], 1);
      varMat.fill(1.0);
      sigma0[k] = newSigma0 * varMat;
      sigma1[k] = newSigma1 * varMat;
      
    }
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Sample sigma from half-Cauchy posterior distribution */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    
    // Compute residual vector
    Res = Y;
    for(j = 0; j < nGroups; j++){
      Res = Res - as<mat>(Xtilde[j]) * as<vec>(getRecentDoubleS(betaTilde[j])); 
    }
    resInnerProd = (trans(Res) * Res);
    // Rcout << "resInnerProd " << resInnerProd(0) << std::endl;
    
    // Sample sigma
    newSigma = sampleSigmaCauchyS(rigamma, n, resInnerProd(0), S, currentSigma);
    //newSigma = rigammaCS(rigamma,  (n + 0.0) / 2.0, (resInnerProd(0) + 1.0) /2.0);
    //newSigma = sqrt(newSigma);
    
    sigmaTildeVec.push_back(newSigma);

    // Save variables if the index is > burn in period
    if(i > BurnIn * outerIter){
      etaCumulative = etaCumulative + etaCount;
      for(j = 0; j < nGroups; j++){
        g1 = gammaCumulative[j];
        g2 = gammaTilde[j];
        gammaCumulative[j] = g1 + g2;
      }
    }
    
  }
  
  // Compute final (output) objects
  validLoops = outerIter - BurnIn * outerIter;
  etaFinal   = as<NumericVector>(etaCumulative) / validLoops;
  for(j = 0; j < nGroups; j++){
    gammaFinal[j] = (as<NumericVector>(gammaCumulative[j])) / validLoops;
  }
  
  return List::create(
    Named("etaFinal") = etaFinal,
    Named("gammaFinal") = gammaFinal,
    Named("betaNew") = betaTilde,
    Named("sigmaTildeVec") = sigmaTildeVec, 
    Named("sigma0") = sigma0, 
    Named("sigma1") = sigma1, 
    Named("theta") = theta
  );

}

