#include <TMB.hpp>
//--------------------------------------------------------
//QGLLVM
//Author: Bert van der Veen
//------------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  //declares all data and parameters used
  DATA_MATRIX(y);
  DATA_MATRIX(x);
  DATA_MATRIX(xr);
  DATA_MATRIX(offset);
  
  //DATA_INTEGER(int_n);//number of quadrature points for Romberg integration in the future
  //DATA_INTEGER(n_int);
  PARAMETER_MATRIX(r0);
  PARAMETER_MATRIX(b);
  PARAMETER_MATRIX(B);
  PARAMETER_VECTOR(lambda);
  PARAMETER_MATRIX(lambda2);
  PARAMETER_VECTOR(lambda3);
  
  //latent variables, u, are treated as parameters
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(lg_phi);
  PARAMETER(log_sigma);// log(SD for row effect)
  DATA_SCALAR(gamma);
  DATA_UPDATE(gamma); //allows updating on R side without recompilation so I can write a function for regularization. DONT IFELSE. Currently doesn't work.
  DATA_SCALAR(gamma2);
  DATA_UPDATE(gamma2);
  DATA_VECTOR(theta4);
  DATA_INTEGER(num_lv);
  DATA_INTEGER(family);
  DATA_INTEGER(max);
  
  PARAMETER_VECTOR(Au);
  PARAMETER_VECTOR(lg_Ar);
  PARAMETER_VECTOR(zeta);
  //DATA_SCALAR(extra);
  //DATA_INTEGER(method);// 0=VA, 1=LA
  DATA_INTEGER(model);
  DATA_VECTOR(random);//random row
  DATA_INTEGER(zetastruc);
  
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
      // if (j == i){
      //   newlam(i,j) = (newlam(i,j));
      // }
    }
  }
  
  //parallel_accumulator<Type> nll(this); // initial value of log-likelihood
  Type nll = 0;
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
  if(max==0){
  //A is a num.lv x num.lv x n array, theta is p x num.lv matrix
  for (int i=0; i<n; i++) {
    nll -= 0.5*(log(A.col(i).matrix().determinant()) - A.col(i).matrix().trace());// log(det(A_i))-sum(trace(A_i))*0.5 sum.diag(A)
  }
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
  
  
  
  matrix <Type> newlam2(num_lv,p);
    for (int j=0; j<lambda2.cols(); j++){
      for (int q=0; q<num_lv; q++){
        newlam2(q,j) = fabs(lambda2(q,j)) + fabs(lambda3(q)) + theta4(q);
      }
    } 

  array<Type> D(num_lv,num_lv,p);
  D.fill(0.0);
  for (int j=0; j<p; j++){
    for (int q=0; q<num_lv; q++){
      D(q,q,j) = 2*newlam2(q,j);
    }
  }
  matrix <Type> eta = C + u*newlam;
  //add exspectation quadratic term
  for (int j=0; j<p; j++){
    for (int i=0; i<n; i++){
      eta(i,j) += -0.5*((u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value());
      if(max==0){
        eta(i,j) += -0.5*((D.col(j).matrix()*A.col(i).matrix()).trace());
      }    
    
    }
  }
  
  if(family==0){
    //likelihood
    if(max==0){
    matrix <Type> e_eta(n,p);
    e_eta.fill(0.0);
    matrix <Type> B(num_lv,num_lv);
    matrix <Type> v(num_lv,1);
    for (int i=0; i<n; i++) {
      matrix <Type> Q = atomic::matinv(A.col(i).matrix());
      for (int j=0; j<p;j++){
        B = (D.col(j).matrix()+Q);
        v = (newlam.col(j)+Q*u.row(i).transpose());
        Type detB = pow(B.determinant(),-0.5);
        Type detA = pow(A.col(i).matrix().determinant(),-0.5);
        e_eta(i,j) += exp(C(i,j) + 0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
        
        nll -= y(i,j)*eta(i,j) - e_eta(i,j) - lfactorial(y(i,j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    }else{
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
       nll -= dpois(y(i,j),exp(eta(i,j)),true);
        }
      }
    }
     }else if(family==1){
       if(max==0){
       matrix <Type> e_eta(n,p);
       e_eta.fill(0.0);
       matrix <Type> B(num_lv,num_lv);
       matrix <Type> v(num_lv,1);
       for (int i=0; i<n; i++) {
         matrix <Type> Q = atomic::matinv(A.col(i).matrix());
         for (int j=0; j<p;j++){
           B = (-D.col(j).matrix()+Q);
           v = (-newlam.col(j)+Q*u.row(i).transpose());
           Type detB = pow(B.determinant(),-0.5);
           Type detA = pow(A.col(i).matrix().determinant(),-0.5);
           e_eta(i,j) += exp((u.row(i)*newlam.col(j)).value()-(u.row(i)*D.col(j).matrix()*u.row(i).transpose()).value()-(D.col(j).matrix()*A.col(i).matrix()).trace()+0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
           
         nll -= y(i,j)*(eta(i,j)-e_eta(i,j)) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j)-e_eta(i,j))) + lgamma(y(i,j)+iphi(j)) - iphi(j)*e_eta(i,j) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
         }
         nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
       }
       }else{
         for (int i=0; i<n; i++) {
           for (int j=0; j<p;j++){
             nll -= dnbinom(y(i,j),1-(exp(eta(i,j))/(1/iphi(j)+exp(eta(i,j)))), 1/iphi(j), true);
           }
         }
         
       }
  } else if(family==2){
    matrix <Type> mu(n,p);
    if(max==0){
    for (int i=0; i<n; i++) {
      for (int j=0; j<p;j++){
        mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
        nll -= dbinom(y(i,j),Type(1),mu(i,j),true);//line below differs from pdf as D is positive only.
        nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - 0.25*(D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 0.5*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() + (u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    }else{
      for (int i=0; i<n; i++) {
        for (int j=0; j<p;j++){
          mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
          nll -= dbinom(y(i,j),Type(1),mu(i,j),true);
        }
      }
    }
  } else if(family==7 && zetastruc==1){
    
    int ymax =  CppAD::Integer(y.maxCoeff());
    int K = ymax - 1;
    
    matrix <Type> zetanew(p,K);
    zetanew.fill(0.0);
    
    int idx = 0;
    for(int j=0; j<p; j++){
      int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
      int Kj = ymaxj - 1;
      if(Kj>1){
        for(int k=0; k<(Kj-1); k++){
          if(k==1){
            zetanew(j,k+1) = fabs(zeta(idx+k));//second cutoffs must be positive
          }else{
            zetanew(j,k+1) = zeta(idx+k);
          }
          
        }
      }
      idx += Kj-1;
    }
    
    for (int i=0; i<n; i++) {
      for(int j=0; j<p; j++){
        int ymaxj = CppAD::Integer(y.col(j).maxCoeff());
        //minimum category
        if(y(i,j)==1){
          nll -= log(pnorm(zetanew(j,0) - eta(i,j), Type(0), Type(1)));
        }else if(y(i,j)==ymaxj){
          //maximum category
          int idx = ymaxj-2;
          nll -= log(1 - pnorm(zetanew(j,idx) - eta(i,j), Type(0), Type(1)));
        }else if(ymaxj>2){
          for (int l=2; l<ymaxj; l++) {
            if(y(i,j)==l && l != ymaxj){
              nll -= log(pnorm(zetanew(j,l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(j,l-2)-eta(i,j), Type(0), Type(1)));
            }
          }
        }
        
        //line below is different as D is positive only and *2
       if(max==0) nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - 0.25*(D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 0.5*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() + (u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
        //nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
      }
      if(max==0) nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    
  } else if(family==7 && zetastruc==0){
    
    int ymax =  CppAD::Integer(y.maxCoeff());
    int K = ymax - 1;
    
    vector <Type> zetanew(K);
    zetanew.fill(0.0);
    for(int k=0; k<(K-1); k++){
      if(k==1){
        zetanew(k+1) = fabs(zeta(k));//second cutoffs must be positive
      }else{
        zetanew(k+1) = zeta(k);
      }
    }
    for (int i=0; i<n; i++) {
      for(int j=0; j<p; j++){
        //minimum category
        if(y(i,j)==1){
          nll -= log(pnorm(zetanew(0) - eta(i,j), Type(0), Type(1)));
        }else if(y(i,j)==ymax){
          //maximum category
          int idx = ymax-2;
          nll -= log(1 - pnorm(zetanew(idx) - eta(i,j), Type(0), Type(1)));
        }else if(ymax>2){
          for (int l=2; l<ymax; l++) {
            if(y(i,j)==l && l != ymax){
              nll -= log(pnorm(zetanew(l-1)-eta(i,j), Type(0), Type(1))-pnorm(zetanew(l-2)-eta(i,j), Type(0), Type(1)));
            }
          }
        }//line below is different as D is positive only and *2
        if(max==0) nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - 0.25*(D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 0.5*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() + (u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
        //  nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
      }
      if(max==0) nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
  }else if(family==3){
    matrix <Type> eta2(n,p);
    if(max==0){
    for (int i=0; i<n; i++) {
      for (int j=0; j<p;j++){
        eta2(i,j) = (newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() + 0.5*(D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() + (u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
        //nll -= -0.5*log(iphi(j)*iphi(j)) - 1/(2*iphi(j)*iphi(j))*(pow(y(i,j)-eta(i,j),2) + eta2(i,j));
        nll -= -0.5*log(iphi(j)*iphi(j)) - 1/(2*iphi(j)*iphi(j))*(y(i,j)*y(i,j) -2*y(i,j)*eta(i,j) + eta(i,j)*eta(i,j) + eta2(i,j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    }else{
      for (int i=0; i<n; i++) {
      for (int j=0; j<p;j++){
      nll -= dnorm(y(i,j), eta(i,j), iphi(j), true);
      }
    }
    }
  }else if(family==4){
    if(max==0){
    matrix <Type> e_eta(n,p);
    e_eta.fill(0.0);
    matrix <Type> B(num_lv,num_lv);
    matrix <Type> v(num_lv,1);
    for (int i=0; i<n; i++) {
      matrix <Type> Q = atomic::matinv(A.col(i).matrix());
      for (int j=0; j<p;j++){
        B = (-D.col(j).matrix()+Q);
        v = (-newlam.col(j)+Q*u.row(i).transpose());
        Type detB = pow(B.determinant(),-0.5);
        Type detA = pow(A.col(i).matrix().determinant(),-0.5);
        e_eta(i,j) += exp(-C(i,j) + 0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
        
        nll -=  ( -eta(i,j) - e_eta(i,j)*y(i,j) )/iphi(j) + log(y(i,j)/iphi(j))/iphi(j) - log(y(i,j)) -lgamma(1/iphi(j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
  }
  }else{
    for (int i=0; i<n; i++) {
      for (int j=0; j<p;j++){
        nll -= dgamma(y(i,j), 1/iphi(j), iphi(j)*exp(eta(i,j)), true);
      }
    }
      
  }

  //shrinks LVs, linear ridge
    for(int q=0; q<(num_lv-1); q++){
      nll += (newlam.row(q).array()*newlam.row(q).array()+newlam2.row(q).array()*newlam2.row(q).array() + newlam.row(q+1).array()*newlam.row(q+1).array()+newlam2.row(q+1).array()*newlam2.row(q+1).array()).sum()*gamma;
      nll += (newlam.row(q+1).array()*newlam.row(q+1).array()+newlam2.row(q+1).array()*newlam2.row(q+1).array()).sum()*gamma;
    }
//shrinks LVs, quadratic ridge
    for(int q=0; q<num_lv; q++){
      nll += (lambda2.row(q).array()*lambda2.row(q).array()).sum()*gamma2;
    }
   if(max==0) nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i

  return nll;
}


