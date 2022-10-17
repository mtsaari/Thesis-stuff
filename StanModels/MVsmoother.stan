functions{
  matrix three_block(matrix m1, matrix m2, matrix m3){
    // creates a block diagonal matrix from  square matrices
    int d1 = rows(m1);
    int d2 = rows(m2);
    int d3 = rows(m3);
    int dmm = d1+d2+d3;
    matrix[dmm,dmm] mm = rep_matrix(0,dmm,dmm);
    mm[1:d1,1:d1] = m1;
    mm[(d1+1):(d1+d2),(d1+1):(d1+d2)] = m2;
    mm[(dmm-d3+1):dmm,(dmm-d3+1):dmm] = m3;
    return mm;
  }
  matrix id_kron_prod(matrix m, int D) {
    int r = rows(m);
    matrix[D*r,D*r] Mprod = rep_matrix(0,D*r,D*r); 
    for(i in 1:D){
      Mprod[(i-1)*r+1:i*r,(i-1)*r+1:i*r] = m;
    }
    return Mprod;
  }
  matrix kron_prod(matrix m1, matrix m2) {
    int r1 = rows(m1);
    int c1 = cols(m1);
    int r2 = rows(m2); 
    int c2 = cols(m2);
    matrix[r1*r2,c1*c2] mprod;
    for(i in 1:r1){
      for(j in 1:c1){
        mprod[((i-1)*r2+1):i*r2, ((j-1)*c2+1):j*c2] = m1[i,j]*m2;
      }
    }
    return mprod;
  }
  matrix Fmatrix(real ell, int D) {
    // creates matrix F corresponding to a Matern kernel with nu = D - 1/2 
    matrix[D,D] F;
    real lambda = sqrt(2*(D-0.5))/ell;
    for (j in 1:D) {
      F[D,j] = -choose(D,j-1)*lambda^(D-j+1);
      for (i in 1:(D-1)) {
        if(j == i+1)
          F[i,j] = 1;
        else
          F[i,j] = 0;
      }
    }
    return F;
  }
  matrix Pmatrix(matrix F){
    // computes steady state covariance matrix corresponding to F
    int D = rows(F);
    real lambda = -F[D,D]/D;
    matrix[D,D] P;
    vector[D*D] Pvec;
    matrix[D,D] ID = diag_matrix(rep_vector(1,D));
    vector[D*D] Q = rep_vector(0,D*D);
    Q[D*D] = -tgamma(D)^2 / tgamma(2*D-1) * (2*lambda)^(2*D-1);
    Pvec = (kron_prod(ID,F)+kron_prod(F,ID))\Q;
    for (j in 1:D) {
      P[1:D,j] = Pvec[(j-1)*D+1:j*D];
    }
    return P;
  }
  matrix Fseasonal(data real per,int J){
    real omega = 2*pi()/per;
    matrix[2*J,2*J] F = rep_matrix(0,2*J,2*J);
    for(j in 1:J){
      F[2*j,2*j-1] = omega*j;
      F[2*j-1,2*j] = -F[2*j,2*j-1];
    }
    return F;
  }
  matrix Pseas(int J, real ell){
    matrix[2*J,2*J] P= rep_matrix(0,2*J,2*J);
    for(i in 1:J){
      int n = (J-i)/2;
      for(j in 0:n){
        P[2*i,2*i] += (2*ell^2)^(-i-2*j)/(tgamma(i+j+1)*tgamma(j+1));
      }
      P[2*i-1,2*i-1] = P[2*i,2*i];
    }
    return 2*exp(-ell^(-2))*P;
  }
}

data{
  int N; // n of time points
  int Nobs; // n of observations
  int D; // n of parallel time series
  int J; // degree of cosine covariance function
  int np; // prediction days
  vector[D] y[N]; // data
  int obs_row[N]; // array for n of observations/time point
  int ri[Nobs]; // array of D-indices 
  
}
transformed data{
  int be[N];
  int sdim = D*(5+2*J);
  matrix[D,D] Im = diag_matrix(rep_vector(1,D));
  matrix[D,sdim] H;
  row_vector[5+J*2] Hvec = rep_row_vector(0,5+J*2);
  matrix[J*2,J*2] Fseas = Fseasonal(7, J);
  Hvec[1:5] = [1,0,1,0,0];
  be[1] = 1;
  for(i in 2:N) {
    be[i] = be[i-1]+obs_row[i-1];
  }
  for(i in 1:J){
    Hvec[4+2*i] = 1;
  }
  H = kron_prod(Im, to_matrix(Hvec,1,5+J*2));
}
parameters{
  vector[D] mu;
  real<lower=0> ellshort;
  real<lower=0> elllong;
  real<lower=0> ellseas;
  real<lower=0> ells2;
  real<lower=0> s2short;
  real<lower=0> s2long;
  vector<lower=0>[D] sspat;
  real<lower=0> sper;
  cholesky_factor_corr[D] L;
  vector<lower=0>[D] epsilon;
}
transformed parameters{
  matrix[2,2] Fshort = Fmatrix(ellshort,2);
  matrix[3,3] Flong = Fmatrix(elllong,3);
  matrix[2*J,2*J] Fper = Fseas - diag_matrix(rep_vector(1/ells2,2*J));
  matrix[D,D] Kc = multiply_lower_tri_self_transpose(diag_pre_multiply(sspat,L));
  matrix[sdim,sdim] Pinf = kron_prod(Kc,three_block(s2short*Pmatrix(Fshort), s2long*Pmatrix(Flong),sper*Pseas(J,ellseas)));
  matrix[sdim,sdim] At = id_kron_prod(matrix_exp(three_block(Fshort,Flong,Fper)),D);
  matrix[sdim,sdim] P = Pinf;
  matrix[sdim,sdim] Q = Pinf - quad_form(Pinf,At');
  matrix[sdim,D] K; matrix[D,D] S; vector[D] v;
  vector[sdim] m = rep_vector(0,sdim); 
  real logl = 0;
  matrix[D,D] meps = diag_matrix(epsilon);
  matrix[sdim,sdim] II = diag_matrix(rep_vector(1,sdim));
  matrix[sdim,sdim] Parr[N];
  matrix[sdim,sdim] Larr[N];
  vector[sdim] marr[N];
  vector[D] invSv = rep_vector(0,D);
  vector[D] varr[N];
    // Kalman filter 
  for (i in 1:N) {
    int id[obs_row[i]] = segment(ri,be[i],obs_row[i]);
    m = At*m; 
    P = quad_form(P,At') + Q;
    Parr[i] = P; // save state cov for smoothing
    marr[i] = m; // save mean for smoothing
    varr[i] = rep_vector(0.0,D);
    if(obs_row[i]>0){
      K = P*H';
      S[id,id] = H[id,:]*K[:,id] + meps[id,id];
      v[id] = y[i][id]-H[id,:]*m - mu[id];
      invSv[id] = S[id,id]\v[id];
      K[:,id] = K[:,id]/S[id,id];
      m += K[:,id]*v[id];
      P -= quad_form(S[id,id],K[:,id]');
      logl += -0.5*(log_determinant(S[id,id]) + dot_product(v[id],invSv[id]));
      }
    Larr[i] = At*(II-(K[:,id]*H[id,:])); // save matrix L for smoothing
    varr[i] = invSv;
  }
} 
model{
  mu ~ normal(0,5);
  ellshort ~ inv_gamma(4,20);
  elllong ~ inv_gamma(5,300);
  ellseas ~ inv_gamma(3,2);
  ells2 ~ inv_gamma(5,300.0);
  L ~ lkj_corr_cholesky(1.5);
  epsilon ~ inv_gamma(3,2);
  s2short ~ inv_gamma(3,2);
  s2long ~ inv_gamma(3,2);
  sspat ~ inv_gamma(3,2);
  sper ~ inv_gamma(3,2);
  target += logl;
}
generated quantities{
 vector[sdim] r = rep_vector(0,sdim);
 vector[sdim] msmooth;
 vector[D] smooth_mean[N];
 vector[D] smooth_short[N];
 vector[D] smooth_long[N];
 vector[D] smooth_per[N];
 vector[D] preds[np];
 vector[D] sims[np];
 matrix[D,D] Covs[np];
 matrix[sdim,sdim] PT = P;
 matrix[sdim,sdim] Pp = P;
 vector[sdim] mp = m;
 for(i in 0:(N-1)){
   int id[obs_row[N-i]] = segment(ri,be[N-i],obs_row[N-i]);
   r = H[id,:]'*varr[N-i][id] + Larr[N-i]'*r;
   msmooth = marr[N-i] + Parr[N-i]*r; 
   smooth_mean[N-i] = H*msmooth + mu;
   for(j in 1:D){
     smooth_short[N-i][j] = msmooth[(j-1)*(sdim/D)+1];
     smooth_long[N-i][j] = msmooth[(j-1)*(sdim/D)+3];
   }
   smooth_per[N-i] = smooth_mean[N-i]-smooth_short[N-i]-smooth_long[N-i]-mu;
 }
 for(k in 1:np){
   mp = At*mp;
   Pp = At*Pp*At' + Q;
   preds[k] = H*mp + mu;
   Covs[k] = H*Pp*H' + meps;
   sims[k] = multi_normal_rng(preds[k],Covs[k]);
 }
}