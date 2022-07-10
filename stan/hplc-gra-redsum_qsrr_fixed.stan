functions {

// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
vector lower_tri(matrix mat) {

int d = rows(mat);
int lower_tri_d = d * (d - 1) / 2;
vector[lower_tri_d] lower;
int count = 1;
for(r in 2:d) {
for(c in 1:(r - 1)) {
lower[count] = mat[r,c];
count += 1;
}
}
return(lower); 
}

// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
real lkj_corr_point_lower_tri_lpdf(matrix rho, vector point_mu_lower, vector point_scale_lower) {

real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(lower_tri(rho) | point_mu_lower, point_scale_lower);
return(lpdf);
}


real lkj_corr_cholesky_point_lower_tri_two_lpdf(matrix cor_L, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 1);
    int d = rows(cor_L);
    matrix[d,d] cor = multiply_lower_tri_self_transpose(cor_L);
    lpdf += normal_lpdf(cor[2,1] | point_mu_lower, point_scale_lower);
    return(lpdf);
  }

// pH and fi at a given time at column inlet
vector gra_state(real t,  vector hplcparam) {

vector[2] sol;
real tg = hplcparam[1];
real td = hplcparam[2];
real fio = hplcparam[5];
real fik = hplcparam[6];
real pHo = hplcparam[8];
real alpha1 = hplcparam[9];
real alpha2 = hplcparam[10];
real fi;

fi = fio+(fik-fio)/tg*(t-td);

if (t<td)
fi = fio;
else if (t>tg+td)
fi = fik;

sol[1]=fi;
sol[2]=pHo+alpha1*fi+alpha2*fi.^2;

return sol;
}

real funlogki(vector logkwx, vector S1, vector pKaw, vector alpha, real S2, vector apH,
              int nDiss, vector chargesA, vector chargesB, vector fipH) {

real logki;
vector[3] logkix;
vector[2] pHmpKa;
real fi=fipH[1];
real pH=fipH[2];

logkix=logkwx-S1*fi/(1+S2*fi)+ chargesA*apH[1]*(pH-7) + chargesB*apH[2]*(pH-7);
pHmpKa=pH-(pKaw+alpha*fi);

if (nDiss==0) {
    logki = logkix[1]; 
}
else if (nDiss==1){
    logki=logkix[1] +
    log1p_exp(log(10)*(pHmpKa[1]+logkix[2]-logkix[1]))/log(10)-
    log1p_exp(log(10)*(pHmpKa[1]))/log(10);
}
else if (nDiss==2){
    logki = logkix[1] +
    log1p_exp(log(10)*(pHmpKa[1]+logkix[2]-logkix[1]) + 
    log1p_exp(log(10)*(pHmpKa[2]+logkix[3]-logkix[2])))/log(10)-
    log1p_exp(log(10)*(pHmpKa[1]) + 
    log1p_exp(log(10)*(pHmpKa[2])))/log(10);
}

return logki;
}

vector areaandslope(real time1, real time2, real invki1, real invki2) {

vector[2] cki_b;
real bo;
real cki;

if (invki2>1.001*invki1) {
    bo = (log(invki2)-log(invki1))/(time2-time1);
    cki = (invki2-invki1)/bo;
}
else {
    bo  = 0.001/(time2-time1);
    cki = (time2-time1)*(invki2+invki1)/2;
}

cki_b[1] = cki;
cki_b[2] = bo;

return cki_b;
}

real chromgratrapz(int steps, 
           vector logkwx, vector logkmx, vector pKaw, vector alpha,
           real S2, vector apH, vector chargesA, vector chargesB,
           int nDiss, vector hplcparam) {

real tg = hplcparam[1];
real td = hplcparam[2];
real to = hplcparam[3];

vector[1] sol;
real time1;
real time2; 
vector[2] fipH1;
vector[2] fipH2;
real logki1;
real logki2; 
real invki1;
real invki2;
vector[2] cki_b;
real cumki1;
real cumki2;
real bo;
real tr;
real dt;

dt = tg/steps;

time1 = 0;
time2 = td;

fipH1 = gra_state(time1,  hplcparam);
// fipH2 = fipH1;

logki1 = funlogki(logkwx, logkmx, pKaw, alpha, S2, apH, nDiss, chargesA, chargesB, fipH1);
//logki2 = logki1;

invki1 = 1/to/10^logki1;
invki2 = invki1;

cumki1 = 0;
cumki2 = td*invki1; // cumulative area

bo     = 0.001/td;  // slope

for(x in 1:steps){ 
    if (cumki2>=1)  continue;
    time1 = time2;
    time2 += dt;
//    fipH1 = fipH2;
    fipH2 = gra_state(time2,  hplcparam);
//    logki1 = logki2;
    logki2 = funlogki(logkwx, logkmx, pKaw, alpha, S2, apH, nDiss, chargesA, chargesB, fipH2);
    invki1 = invki2;
    invki2 = 1/to/10^logki2;
    cki_b = areaandslope(time1, time2, invki1, invki2);
    cumki1 = cumki2;
    cumki2 += cki_b[1]; // cumulative area
    bo      = cki_b[2]; //slope
}

if (cumki2>=1) {
    tr = time1+log1p((1-cumki1)*bo/invki1)/bo;
}
else if (cumki2<1) {
    tr = time2+(1-cumki2)/invki2;
}

return tr;
}

real partial_sum(int[] ind, int start, int end, vector trObs, 
int[] mod, int[] steps, int[] analyte, int[] pHid,
vector[] hplcparam, vector[] chargesA, vector[] chargesB, int[] nDiss,
vector[] logkwx, vector dlogkT, vector[] S1mx, vector[] S1ax,
vector[] pKaw, vector[] alpham, vector[] alphaa,
real S2m, real S2a,
vector apH,
vector sigma) {

real lp = 0;
real y_hat;

for(z in start:end){

if (mod[z]==1) {
y_hat = chromgratrapz(steps[z], 
       logkwx[analyte[z],] + dlogkT[analyte[z]]*hplcparam[z,11], 
       S1mx[analyte[z],],  
       pKaw[analyte[z],], 
       alpham[analyte[z],],
       S2m,
       apH,
       chargesA[analyte[z],], 
       chargesB[analyte[z],], 
       nDiss[analyte[z]],
       hplcparam[z]);
 }

if (mod[z]==2) {
y_hat = chromgratrapz(steps[z], 
       logkwx[analyte[z],]  + dlogkT[analyte[z]]*hplcparam[z,11], 
       S1ax[analyte[z],],  
       pKaw[analyte[z],], 
       alphaa[analyte[z],],
       S2a,
       apH,
       chargesA[analyte[z],],
       chargesB[analyte[z],], 
       nDiss[analyte[z]],
       hplcparam[z]);
}

  real trHat = hplcparam[z,3] + hplcparam[z,4] + y_hat; 

  lp = lp + student_t_lpdf(trObs[z] | 3, trHat, sigma[analyte[z]]);

 }
 return lp;
}
}

data{
int nAnalytes;	           // number of analytes
int nObs;		           // number of observations
int npH;                   // npH;
int analyte[nObs];	       // analyte indexes
int pHid[nObs];
int<lower=1> steps[nObs];   // steps for gradient retention time aproimation
vector[11] hplcparam[nObs]; // [tg, td, to, te, fio, fik, mod, pHo, alpha1, alpha2, (temp-25)/10]
int<lower=0> mod[nObs];     // MeOH==1, ACN==2 (repeats hplcparam(:,7))

vector[nAnalytes] logPobs; 

int<lower=0,upper=2> maxR;
int<lower=0,upper=2> R[nAnalytes];
ordered[maxR] pKaslit[nAnalytes];
vector[maxR] pKasliterror[nAnalytes];
vector[maxR] groupsA[nAnalytes];
vector[maxR] groupsB[nAnalytes];
vector[maxR+1] chargesA[nAnalytes];
vector[maxR+1] chargesB[nAnalytes];

int<lower=0> K;                      //  number of predictors (functional groups)
matrix[nAnalytes, K] nrfungroups;   // predictor matrix (functional groups)   

vector[nObs] trobs; // observed retention factors 

vector[K] mpilogkw;
vector[K] spilogkw;
vector[K] mpiS1m;
vector[K] spiS1m;
vector[K] mpiS1a;
vector[K] spiS1a;
}

transformed data {
int grainsize = 1;
int ind[nObs] = rep_array(1, nObs);
vector[3] point_mu_lower = [0.7811,0.7135,0.9148]';       // mean priors for rho
vector[3] point_scale_lower = [0.0395,0.0442,0.0176]'; // std priors for rho
}

parameters{
real logkwHat;	       // typical value of logkw [N]
real S1mHat;	       // typical value of S1m [N]
real S1aHat;           // typical value of S1a [N]
real dlogkwHat[2];     // typical value of dlogkw [A,B] 
real dSmHat[2];        // typical value of dlogkm [A,B] 
real dSaHat[2];        // typical value of dlogka [A,B] 
real<lower = 0> S2mHat; // typical value of S2m 
real<lower = 0> S2aHat; // typical value of S2a
vector[3] beta;         // effects of logP 
real dlogkTHat;         // typical dlogkT
vector[2] alphaAHat;  // changes of pKa with org. mod for acids [MeOH, ACN]
vector[2] alphaBHat;  // changes of pKa with org. mod for bases [MeOH, ACN]

vector<lower = 0.01>[3] omega;   // between analyte variabilities (neutral forms)
corr_matrix[3] rho1;	                        // correlation matrix	 
vector<lower = 0.01>[3] kappa;    // between analyte variabilities (diss. forms)
vector<lower = 0.01>[2] tau;     // between analyte variabilities for acids pKa
cholesky_factor_corr[2] L2;	                        // cholesky 
real<lower = 0.01, upper = 1> omegadlogkT;  // between analyte variability for temperature

// between buffer differences
vector[2] apH; // pH effects

vector[K] pilogkw;  // regression coefficient for logkw
vector[K] piS1m;  // regression coefficient for S1m
vector[K] piS1a;  // regression coefficient for S1a

vector<lower = 0.01>[3] sdpi;     // between analyte variabilities for acids pKa

// residual variability
real<lower = 0.01> msigma; // mean
real<lower = 0.01> ssigma; // scale

// individual values of chromatographic parameters
vector[3] param[nAnalytes]; 
vector[nAnalytes] dlogkT;	
matrix[nAnalytes,maxR+1] dlogkwA;
matrix[nAnalytes,maxR+1] dlogkwB;
matrix[nAnalytes,maxR+1] dSmA;
matrix[nAnalytes,maxR+1] dSmB;
matrix[nAnalytes,maxR+1] dSaA;
matrix[nAnalytes,maxR+1] dSaB;

vector[maxR] pKaw[nAnalytes];
matrix[maxR,nAnalytes] etaStd1;
matrix[maxR,nAnalytes] etaStd2;

// and residuals
vector<lower = 0.01, upper = 4>[nAnalytes] sigma;
}

transformed parameters{
vector[maxR+1] logkwx[nAnalytes];
vector[maxR+1] S1mx[nAnalytes];
vector[maxR+1] S1ax[nAnalytes];
matrix[nAnalytes,maxR] alpha1;   //MeOH or ACN
matrix[nAnalytes,maxR] alpha2;   //MeOH or ACN

vector[maxR] alpham[nAnalytes];
vector[maxR] alphaa[nAnalytes];

vector[3] miu[nAnalytes];	
cov_matrix[3] Omega; // variance-covariance matrix

Omega = quad_form_diag(rho1, omega);	// diag_matrix(omega) * rho * diag_matrix(omega)


// Matt's trick to use unit scale 
 alpha1 = diag_pre_multiply(tau, L2 * etaStd1)';
 alpha2 = diag_pre_multiply(tau, L2 * etaStd2)';


for(i in 1:nAnalytes){
    miu[i,1]  = logkwHat + beta[1] * (logPobs[i]-2.2) + nrfungroups[i,1:K] * pilogkw;
    miu[i,2]  = S1mHat   + beta[2] * (logPobs[i]-2.2) + nrfungroups[i,1:K] * piS1m; 
    miu[i,3]  = S1aHat   + beta[3] * (logPobs[i]-2.2) + nrfungroups[i,1:K] * piS1a;
}

for(i in 1:nAnalytes){
for(r in 1:maxR+1){
logkwx[i,r] = param[i, 1] +
            dlogkwA[i,r]*chargesA[i,r] +
            dlogkwB[i,r]*chargesB[i,r];
S1mx[i,r] = (param[i, 2] + 
            dSmA[i,r]*chargesA[i,r] +
            dSmB[i,r]*chargesB[i,r])*(1+S2mHat);
S1ax[i,r] = (param[i, 3] + 
            dSaA[i,r]*chargesA[i,r] +
            dSaB[i,r]*chargesB[i,r])*(1+S2aHat);

}}

for(i in 1:nAnalytes){
alpham[i,1] = (alphaAHat[1]+alpha1[i,1]) * groupsA[i,1] + (alphaBHat[1]+alpha1[i,1]) * groupsB[i,1];
alpham[i,2] = (alphaAHat[1]+alpha2[i,1]) * groupsA[i,2] + (alphaBHat[1]+alpha2[i,1]) * groupsB[i,2];

alphaa[i,1] = (alphaAHat[2]+alpha1[i,2]) * groupsA[i,1] + (alphaBHat[2]+alpha1[i,2]) * groupsB[i,1];
alphaa[i,2] = (alphaAHat[2]+alpha2[i,2]) * groupsA[i,2] + (alphaBHat[2]+alpha2[i,2]) * groupsB[i,2];
}

}
model{
logkwHat  ~ normal(3.1543,0.0871);
S1mHat    ~ normal(4.5141,0.0971);
S1aHat    ~ normal(5.5600,0.1348);
dlogkwHat[1] ~ normal(-0.7357,0.0588);
dlogkwHat[2] ~ normal(-0.9359,0.0461);
dSmHat[1]    ~ normal(0.3311,0.1030);
dSmHat[2]    ~ normal(0.1098,0.0704);
dSaHat[1]    ~ normal(0.8910,0.1057);
dSaHat[2]    ~ normal(-0.4577,0.0670);
S2mHat    ~ normal(0.3741,0.0250);
S2aHat    ~ normal(0.8194,0.0342);
alphaAHat[1] ~ normal(1.9736,0.1703);
alphaAHat[2] ~ normal(2.1454,0.1844);
alphaBHat[1] ~ normal(-1.0005,0.1405);
alphaBHat[2] ~ normal(-0.8940,0.1641);
beta[1] ~ normal(0.7442,0.0313);
beta[2] ~ normal(0.3553,0.0393);
beta[3] ~ normal(0.3895,0.0513);
omega[1]       ~ normal(0.6150,0.0393);
omega[2]       ~ normal(0.6762,0.0469);
omega[3]       ~ normal(0.9206,0.0631);
rho1         ~ lkj_corr_point_lower_tri(point_mu_lower, point_scale_lower);
kappa[1]       ~ normal(0.5305,0.0276);
kappa[2]       ~ normal(0.5526,0.0447);
kappa[3]       ~ normal(0.5409,0.0422);

apH[1] ~ normal(-0.0238,0.0010);
apH[2] ~ normal( 0.0851,0.0009);

for(k in  1:K){
pilogkw[k] ~ normal(mpilogkw[k],spilogkw[k]);
piS1m[k]   ~ normal(mpiS1m[k],spiS1m[k]);
piS1a[k]   ~ normal(mpiS1a[k],spiS1a[k]);
}

sdpi[1] ~ normal(0.1964,0.0301);
sdpi[2] ~ normal(0.1736,0.0288);
sdpi[3] ~ normal(0.3162, 0.0387);

tau[1] ~ normal(2.2561,0.1644);
tau[2] ~ normal(2.5580,0.1841);

L2 ~ lkj_corr_cholesky_point_lower_tri_two(0.9408, 0.0165);

dlogkTHat   ~ normal(-0.0946,0.0026);
omegadlogkT ~ normal(0.0344,0.0021);

sigma  ~ lognormal(log(msigma),ssigma); 
msigma ~ normal(0.3671,0.0278);
ssigma ~ normal(0.9989,0.0536);

for(i in  1:nAnalytes){
param[i] ~ multi_normal(miu[i],Omega);
}

to_vector(dlogkwA) ~ normal(dlogkwHat[1],kappa[1]);
to_vector(dlogkwB) ~ normal(dlogkwHat[2],kappa[1]);
to_vector(dSmA) ~ normal(dSmHat[1],kappa[2]);
to_vector(dSmB) ~ normal(dSmHat[2],kappa[2]);
to_vector(dSaA) ~ normal(dSaHat[1],kappa[3]);
to_vector(dSaB) ~ normal(dSaHat[2],kappa[3]);

to_vector(etaStd1) ~ std_normal();
to_vector(etaStd2) ~ std_normal();

dlogkT  ~ normal(dlogkTHat,omegadlogkT);

for (i in 1:nAnalytes) pKaw[i] ~ normal(pKaslit[i],pKasliterror[i]);

target += reduce_sum(partial_sum, ind, grainsize, trobs, 
        mod, steps, analyte, pHid, hplcparam, chargesA, chargesB, R,
        logkwx, dlogkT, S1mx, S1ax,
        pKaw, alpham, alphaa,
        S2mHat, S2aHat,
        apH, sigma);
}

generated quantities{
real trCond[nObs];
real trPred[nObs];
real y_hat_Cond;
real y_hat_Pred;
vector[3] paramPred[nAnalytes]; 
matrix[nAnalytes,maxR+1] dlogkwAPred;
matrix[nAnalytes,maxR+1] dlogkwBPred;
matrix[nAnalytes,maxR+1] dSmAPred;
matrix[nAnalytes,maxR+1] dSmBPred;
matrix[nAnalytes,maxR+1] dSaAPred;
matrix[nAnalytes,maxR+1] dSaBPred;
vector[nAnalytes] dlogkTPred;
vector[maxR+1] logkwxPred[nAnalytes];
vector[maxR+1] S1mxPred[nAnalytes];
vector[maxR+1] S1axPred[nAnalytes];
vector[maxR] pKawPred[nAnalytes];
vector[nAnalytes] sigmaPred; 
matrix[maxR,nAnalytes]  etaStd1Pred;
matrix[maxR,nAnalytes]  etaStd2Pred;
matrix[nAnalytes,maxR] alpha1Pred; 
matrix[nAnalytes,maxR] alpha2Pred; 
vector[maxR] alphamPred[nAnalytes];
vector[maxR] alphaaPred[nAnalytes];

corr_matrix[2] rho2;
 
rho2 = L2 * L2';
  
for(i in 1:nAnalytes){
dlogkTPred[i] = normal_rng(dlogkTHat,omegadlogkT); 
sigmaPred[i] = lognormal_rng(log(msigma), ssigma);
}

for(i in 1:nAnalytes){
paramPred[i] = multi_normal_rng(miu[i],Omega);
}

for(r in 1:(maxR+1)){ 
for(i in 1:nAnalytes){
dlogkwAPred[i, r] = normal_rng(dlogkwHat[1],kappa[1]);
dlogkwBPred[i, r] = normal_rng(dlogkwHat[2],kappa[1]);
dSmAPred[i, r] = normal_rng(dSmHat[1],kappa[2]);
dSmBPred[i, r] = normal_rng(dSmHat[2],kappa[2]);
dSaAPred[i, r] = normal_rng(dSaHat[1],kappa[3]);
dSaBPred[i, r] = normal_rng(dSaHat[2],kappa[3]);
}
}

for(r in 1:(maxR)){ 
for(i in 1:nAnalytes){
pKawPred[i,r] = normal_rng(pKaslit[i,r], pKasliterror[i,r]);
etaStd1Pred[r,i] = normal_rng(0,1);
etaStd2Pred[r,i] = normal_rng(0,1);
}
}

alpha1Pred = diag_pre_multiply(tau, L2 * etaStd1Pred)';
alpha2Pred = diag_pre_multiply(tau, L2 * etaStd2Pred)';
 
for(i in 1:nAnalytes){
for(r in 1:maxR+1){
logkwxPred[i,r] = paramPred[i, 1] +
    dlogkwAPred[i,r]*chargesA[i,r] +
    dlogkwBPred[i,r]*chargesB[i,r];
S1mxPred[i,r] = (paramPred[i, 2] + 
    dSmAPred[i,r]*chargesA[i,r] +
    dSmBPred[i,r]*chargesB[i,r])*(1+S2mHat);
S1axPred[i,r] = (paramPred[i, 3] + 
    dSaAPred[i,r]*chargesA[i,r] +
    dSaBPred[i,r]*chargesB[i,r])*(1+S2aHat);
}}

for(i in 1:nAnalytes){
alphamPred[i,1] = (alphaAHat[1]+alpha1Pred[i,1]) * groupsA[i,1] + (alphaBHat[1]+alpha1Pred[i,1]) * groupsB[i,1];
alphamPred[i,2] = (alphaAHat[1]+alpha2Pred[i,1]) * groupsA[i,2] + (alphaBHat[1]+alpha2Pred[i,1]) * groupsB[i,2];

alphaaPred[i,1] = (alphaAHat[2]+alpha1Pred[i,2]) * groupsA[i,1] + (alphaBHat[2]+alpha1Pred[i,2]) * groupsB[i,1];
alphaaPred[i,2] = (alphaAHat[2]+alpha2Pred[i,2]) * groupsA[i,2] + (alphaBHat[2]+alpha2Pred[i,2]) * groupsB[i,2];
}

// COND
 for(z in 1:nObs){

if (mod[z]==1) {
y_hat_Cond = chromgratrapz(steps[z], 
    logkwx[analyte[z],] +  dlogkT[analyte[z]]*hplcparam[z,11], 
    S1mx[analyte[z],],  
    pKaw[analyte[z],], 
    alpham[analyte[z],],
    S2mHat,
    apH,
    chargesA[analyte[z],],
    chargesB[analyte[z],],
    R[analyte[z]],
    hplcparam[z]);
}

if (mod[z]==2) {
y_hat_Cond = chromgratrapz(steps[z], 
    logkwx[analyte[z],] + dlogkT[analyte[z]]*hplcparam[z,11], 
    S1ax[analyte[z],],  
    pKaw[analyte[z],], 
    alphaa[analyte[z],],
    S2aHat,
    apH,
    chargesA[analyte[z],],
    chargesB[analyte[z],],
    R[analyte[z]],
    hplcparam[z]);
}

real trHatCond = hplcparam[z,3] + hplcparam[z,4] + y_hat_Cond; 
trCond[z] = student_t_rng(3,trHatCond, sigma[analyte[z]]);
}

// PRED
for(z in 1:nObs){
if (mod[z]==1) {
y_hat_Pred = chromgratrapz(steps[z], 
       logkwxPred[analyte[z],] +  dlogkTPred[analyte[z]]*hplcparam[z,11], 
       S1mxPred[analyte[z],],  
       pKawPred[analyte[z],], 
       alphamPred[analyte[z],],
       S2mHat,
       apH,
       chargesA[analyte[z],],
       chargesB[analyte[z],],
       R[analyte[z]],
       hplcparam[z]);
}

if (mod[z]==2) {
y_hat_Pred = chromgratrapz(steps[z], 
       logkwxPred[analyte[z],]  + dlogkTPred[analyte[z]]*hplcparam[z,11], 
       S1axPred[analyte[z],],  
       pKawPred[analyte[z],], 
       alphaaPred[analyte[z],],
       S2aHat,
       apH,
       chargesA[analyte[z],],
       chargesB[analyte[z],],
       R[analyte[z]],
       hplcparam[z]);
}

  real trHatPred = hplcparam[z,3] + hplcparam[z,4] + y_hat_Pred; 

  trPred[z] = student_t_rng(3, trHatPred, sigmaPred[analyte[z]]);
}
}
