#define TMB_LIB_INIT R_init_gllvm
#include <TMB.hpp>
#include<math.h>
//--------------------------------------------------------
//GLLVM
//Author: Jenni Niku, Bert van der Veen
//------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  //declares all data and parameters used
  DATA_MATRIX(y);
  DATA_MATRIX(x);
  DATA_MATRIX(xr);
  DATA_MATRIX(offset);
  
  PARAMETER_MATRIX(r0);
  PARAMETER_MATRIX(b);
  PARAMETER_MATRIX(B);
  PARAMETER_VECTOR(lambda);
  PARAMETER_MATRIX(lambda2);
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phi);
  PARAMETER(log_sigma);// log(SD for row effect)
  
  DATA_INTEGER(num_lv);
  DATA_INTEGER(family);
  
  PARAMETER_VECTOR(Au);
  PARAMETER_VECTOR(lg_Ar);
  PARAMETER_MATRIX(zeta);
  //DATA_SCALAR(extra);
  //DATA_INTEGER(method);// 0=VA, 1=LA
  DATA_INTEGER(model);
  DATA_VECTOR(random);//random row
  
  int n = y.rows();
  int p = y.cols();
  
  vector<Type> iphi = exp(lg_phi);
  vector<Type> Ar = exp(lg_Ar);
  Type sigma = exp(log_sigma);
  
  if(random(0)<1){  r0(0,0) = 0;}
  
  matrix<Type> C(n,p);
  C.fill(0.0);
  
  matrix<Type> newlam(num_lv,p);
  
  
  for (int j=0; j<p; j++){
    for (int i=0; i<num_lv; i++){
      if (j < i){
        newlam(i,j) = 0;
      } else{
        newlam(i,j) = lambda(j);
        if (i > 0){
          newlam(i,j) = lambda(i+j+i*p-(i*(i-1))/2-2*i);
        }
      }
      // set diag>0 !!!!!!!!!!!
      if (j == i){
        newlam(i,j) = fabs(newlam(i,j));
      }
    }
  }
  
  Type nll = 0.0; // initial value of log-likelihood
  
  C += r0*xr + offset;
  
  if(random(0)>0){
    for (int i=0; i<n; i++) {
      Ar(i)=pow(Ar(i),2);
    }}
  
  
  array<Type> A(num_lv,num_lv,n);
  for (int d=0; d<(num_lv); d++){
    for(int i=0; i<n; i++){
      A(d,d,i)=exp(Au(d*n+i));
    }
  }
  if(Au.size()>num_lv*n){
    int k=0;
    for (int c=0; c<(num_lv); c++){
      for (int r=c+1; r<(num_lv); r++){
        for(int i=0; i<n; i++){
          A(r,c,i)=Au(num_lv*n+k*n+i);
          A(c,r,i)=A(r,c,i);
        }
        k++;
      }}
  }
  /*Calculates the commonly used (1/2) theta'_j A_i theta_j
  A is a num.lv x nmu.lv x n array, theta is p x num.lv matrix*/
  for (int i=0; i<n; i++) {
    nll -= 0.5*(log(A.col(i).matrix().determinant()) - A.col(i).matrix().trace());// log(det(A_i))-sum(trace(A_i))*0.5 sum.diag(A)
  }
  
  
  
  if(model<1){
    C += x*b;
  } else {
    matrix<Type> eta1=x*B;
    int m=0;
    for (int j=0; j<p;j++){
      for (int i=0; i<n; i++) {
        C(i,j)+=b(0,j)+eta1(m,0);
        m++;
      }
    }
  }
  //lambda2 = lambda2.cwiseAbs(); //sign constraint quadratic effect
  
  matrix <Type> newlam2 = lambda2.cwiseAbs(); //positive only
  matrix <Type> eta = C + u*newlam - (u.array()*u.array()).matrix()*newlam2; //intercept(s), linear effect and negative only quadratic term
  
  array<Type> D(num_lv,num_lv,p);
  D.fill(0.0);
  for (int j=0; j<p; j++){
    for (int q1=0; q1<num_lv; q1++){
      for (int q2=0; q2<num_lv; q2++){
        if(q1!=q2){
          //D(q1,q2,j) = 0.0;
        }else{
          D(q1,q2,j) = 2*newlam2(q1,j);  
        }
      }
    }
  }
  //trace of quadratic effect
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<p;j++){
      eta(i,j) -= 0.5*((D.col(j).matrix()*A.col(i).matrix()).diagonal().sum());//could just remove .matrix().diagonal() and sum all, off-diagonal is zero due to D being a diagonal matrix
    }
  }     
  if(family==0){
    //likelihood
    matrix <Type> e_eta(n,p);
    matrix <Type> B(num_lv,num_lv);
    matrix <Type> v(num_lv,1);
    for (int i=0; i<n; i++) {
      matrix <Type> Q = A.col(i).matrix().inverse();
      for (int j=0; j<p;j++){
        //try to write this like a series of normal multiplication with arrays, -=, += and exp at the end?
        //the problem is in e_eta since y*eta-exp(eta) works fine.
        //its the determinants.
        B = (D.col(j).matrix()+Q);
        v = (newlam.col(j)+Q*u.row(i).transpose());
        Type detB = pow(B.determinant(),-0.5);
        Type detA = pow(A.col(i).matrix().determinant(),-0.5);
        e_eta(i,j) = exp(C(i,j) + 0.5*((v.transpose()*B.inverse()*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
        nll -= y(i,j)*eta(i,j) - e_eta(i,j) - lfactorial(y(i,j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
  }else if(family==1){
    matrix <Type> zetanew(n,p);
    matrix <Type> B(num_lv,num_lv);
    matrix <Type> v(num_lv,1);
    for (int i=0; i<n; i++) {
      matrix <Type> Q = A.col(i).matrix().inverse();
      for (int j=0; j<p;j++){
        B = (D.col(j).matrix()+Q);
        v = (newlam.col(j)+Q*u.row(i).transpose());
        Type detB = pow(B.determinant(),-0.5);
        Type detA = pow(A.col(i).matrix().determinant(),-0.5);
        zetanew(i,j) = iphi(j) + exp(C(i,j) + 0.5*((v.transpose()*B.inverse()*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
        
        nll -= y(i,j) * eta(i,j) - (y(i,j) + iphi(j))*log(zetanew(i,j)) - iphi(j)*((y(i,j) + iphi(j))/zetanew(i,j)) + lgamma(y(i,j)+iphi(j)) - lfactorial(y(i,j)) + iphi(j)*log(iphi(j)) - lgamma(iphi(j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
  } else if(family==2){
    matrix <Type> mu(n,p);
    for (int i=0; i<n; i++) {
      for (int j=0; j<p;j++){
        mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
        //y(i,j)*log(mu(i,j))+(1-y(i,j))*log(1-mu(i,j));
        nll -= dbinom(y(i,j),Type(1),mu(i,j),true);
        nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();   
        //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));// 
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    
  } else if(family==3){
    int ymax = 1;
    int ymin = 999;
    
    //find maximum and minimum integer values in the data
    for (int i=0; i<n; i++) {
      for(int j=0; j<p; j++){
        if(ymax<y(i,j)){
          ymax += 1;
        }
        if(ymax<ymin){
          ymin = ymax;
        }else if(ymin>y(i,j) && ymin != 0){
          ymin -= 1;
        }
      }
    }
    
    //add a column of zeros
    matrix <Type> zetanew(p,zeta.cols()+1);
    zetanew.fill(0.0);
    int K = zeta.cols();
    for (int j=0; j<p; j++){
      for(int k=0; k<K; k++){
        zetanew(j,k+1) = zeta(j,k);
      }
      
    }
    
    for (int i=0; i<n; i++) {
      for(int j=0; j<p; j++){
        //minimum category
        if(y(i,j)==ymin){
          nll -= log(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1)));
        }
        //maxmimum category
        if(y(i,j)==ymax){
          int idx = ymax-(1+ymin);
          nll -= log(1 - pnorm(zetanew(j,idx) - eta(i,j), Type(0), Type(1)));
        }
        //everything else
        for (int l=ymin; l<ymax; l++) {
          if(y(i,j)==l && y(i,j) != ymin && y(i,j) != ymax){
            int idx = (l-ymin);
            int idx2 = l-(1+ymin);
            nll -= log(pnorm(zetanew(j,idx)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(j,idx2)-eta(i,j), Type(0), Type(1))); 
          }
        }
        
        nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();   
        //log(pow(mu(i,j),y(i,j))*pow(1-mu(i,j),(1-y(i,j))));// 
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    
  }
  nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i
  
  
  
  return nll;
}
