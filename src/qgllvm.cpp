#include <TMB.hpp>
//--------------------------------------------------------
//QGLLVM
//Author: Bert van der Veen
//------------------------------------------------------------

template<class Float>
struct integrand {
  typedef Float Scalar; // Required by integrate
  Float C;
  matrix<Float> theta;
  matrix<Float> theta2;
  matrix<Float> u;
  matrix<Float> A;
  matrix<Float> Linv;
  //vector<Float> z;//this doesn't work, breaks the whole thing. The integration variables need to be defined as separate floats.
  Float z1;
  Float z2;
  
  Float operator() () {
    //transform znew to z so we integrate over a std normal and can define the integration limits.
    // Z ~ N(u2,A), znew~N(0,I)
    //Define L as lower cholesky of A^-1
    int num_lv = A.cols();
    matrix <Float> z(1,num_lv);
    z(0,0) = z1;
    z(0,1) = z2;
    //    z << z1, z2;
    
    //vector<Float> z(2);
    
    
    //Transform std normal vector to vector with mean and variance from VA using inverse of cholesky. Not very computationally efficient..
    matrix <Float> znew = u + z*Linv;//.transpose();//no transpose as eigen gives the upper cholesky not the lower
    //need to populate a vector again..
    //vector <Float> u()
    
    
    Float ans = -C;
    //for (int q=0; q<num_lv; q++){
    //  ans -= C + znew(0,q)*theta(0,q) - znew(0,q)*znew(0,q)*theta2(0,q);
    //}
    ans -= (znew.array()*theta.array()).sum() - (znew.array()*znew.array()*theta2.array()).sum();
    //density needs to be evaluated at a vector
    ans -= density::MVNORM(A)(znew.row(0)-u.row(0));
    ans = exp(ans);
    // Float ans = exp(- z(0)*theta(0) - z(1)*theta(1) + z(0)*z(0)*theta2(0)+ z(1)*z(1)*theta2(1)); 
    // ans *= exp(-density::MVNORM(A)(z-u2));
    // Avoid NaNs in the gradient:
  //  if (ans == 0) ans = 0;
    // Avoid NaNs in the tail of the distribution:
  //  using atomic::tiny_ad::isfinite;
  //  if (!isfinite(ans)) ans = 0;
    return ans;
  }
  
  // Integrate wrt (x,y)
  Float integrate() {
    Float ans = gauss_kronrod::mvIntegrate(*this).
    //can write a simple ifelse statement here for # of lVs as wrt doesn't accept a vector
    wrt(z1, -INFINITY, INFINITY).
    wrt(z2, -INFINITY, INFINITY) ();
    return ans;
  }
  
  
};

//wrapper function that takes vector arguments. Then, based on num.lv I pass to various integrate functions
// template<class Float>
// Float eval(vector<Float> input) {
//   int num_lv = CppAD::Integer(input[0]);
//   Float C = input[1];
//   matrix <Float> lambda(0,num_lv);
//   matrix <Float> lambda2(0,num_lv);
//   matrix <Float> u(0,num_lv);
//   matrix <Float> A(num_lv,num_lv);
// //these are matrices because I can pass matrices, but not vectors! Can create vectors though..weird
//     for (int q=0; q<num_lv;q++){
//        lambda(0,q) = input(2+q);
//         lambda2(0,q) = input(2+num_lv+q);
//         u(0,q) = input(2+num_lv*2+q);
//       for (int q2=0; q2<num_lv;q2++){
//         A(q2,q) = input(2+num_lv*3+num_lv*q+q2);
//         //A.col(i).col(q)(q2);
//       }
//     }
//   
//   integrand<Float> f = {num_lv, C, lambda, lambda2, u, A};
//   return f.integrate();
// }//still register atomic if possible, type conversion still a problem for compilation
template<class Float>
Float eval(Float C, matrix<Float> lambda, matrix<Float> lambda2, matrix<Float> u, matrix<Float> A, matrix<Float> Linv) {
  integrand<Float> f = {C, lambda, lambda2, u, A, Linv};
  return f.integrate();
}//might want to register atomic this to reduce the tape size.

template<class Float>
struct integrand2 {
  typedef Float Scalar; // Required by integrate
  matrix<Float> theta;
  matrix<Float> theta2;
  matrix<Float> u;
  matrix<Float> A;
  // matrix<Float> Linv;
  //vector<Float> z;//this doesn't work, breaks the whole thing. The integration variables need to be defined as separate floats.
  Float z1;
  Float z2;
  
  Float operator() () {
    //transform znew to z so we integrate over a std normal and can define the integration limits.
    // Z ~ N(u2,A), znew~N(0,I)
    //Define L as lower cholesky of A^-1
    int num_lv = A.cols();
    matrix <Float> z(1,num_lv);
    z(0,0) = z1;
    z(0,1) = z2;
    //    z << z1, z2;
    
    //vector<Float> z(2);
    
    
    //Transform std normal vector to vector with mean and variance from VA using inverse of cholesky. Not very computationally efficient..
    // matrix <Float> znew = u + z*Linv;//.transpose();//no transpose as eigen gives the upper cholesky not the lower
    //need to populate a vector again..
    //vector <Float> u()
    
    
    Float ans = 0;
    //for (int q=0; q<num_lv; q++){
    //  ans -= C + znew(0,q)*theta(0,q) - znew(0,q)*znew(0,q)*theta2(0,q);
    //}
    ans += (z.array()*theta.array()).sum() - (z.array()*z.array()*theta2.array()).sum();
    //density needs to be evaluated at a vector
    ans -= density::MVNORM(A)(z.row(0)-u.row(0));//could also evaluate this at z with std normal..
    ans = exp(ans);
    // Float ans = exp(- z(0)*theta(0) - z(1)*theta(1) + z(0)*z(0)*theta2(0)+ z(1)*z(1)*theta2(1)); 
    // ans *= exp(-density::MVNORM(A)(z-u2));
    // Avoid NaNs in the gradient:
    //if (ans == 0) ans = 0;
    // Avoid NaNs in the tail of the distribution:
   // using atomic::tiny_ad::isfinite;
    //if (!isfinite(ans)) ans = 0;
    return ans;
  }
  
  // Integrate wrt (x,y)
  Float integrate() {
    using namespace gauss_kronrod;
    // control<Float> c = {100, 1e-6, 1e-6};
    Float ans = mvIntegrate(*this, control(100,1e-4,1e-4)).//100 subdivisions is default, 1e6 is more accurate than default(1e-4)
    //can write a simple ifelse statement here for # of lVs as wrt doesn't accept a vector
    wrt(z1, -INFINITY, INFINITY).
    wrt(z2, -INFINITY, INFINITY) ();
    return ans;
  }
  
  
};

template<class Float>
Float eval2(matrix<Float> lambda, matrix<Float> lambda2, matrix<Float> u, matrix<Float> A){//, matrix<Float> Linv) {
  integrand2<Float> f = {lambda, lambda2, u, A};//, Linv};
  return log(f.integrate());
}//might want to register atomic this to reduce the tape size.


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
  DATA_VECTOR(gamma);
  DATA_MATRIX(gamma2);
  DATA_VECTOR(theta4);
  DATA_INTEGER(num_lv);
  DATA_INTEGER(family);
  
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
  
  
  
  matrix <Type> newlam2(num_lv,p);
  int nTol = lambda2.cols();
  if(nTol==1){
    for (int j=0; j<p; j++){
      for (int q=0; q<num_lv; q++){
        newlam2(q,j) = fabs(lambda3(q)) + theta4(q);
      }
    }  
  }else{
    for (int j=0; j<p; j++){
      for (int q=0; q<num_lv; q++){
        newlam2(q,j) = fabs(lambda2(q,j)) + fabs(lambda3(q)) + theta4(q);
      }
    }
  }
  
  matrix <Type> eta = C + u*newlam - (u.array()*u.array()).matrix()*newlam2; //intercept(s), linear effect and negative only quadratic term
  
  array<Type> D(num_lv,num_lv,p);
  D.fill(0.0);
  for (int j=0; j<p; j++){
    for (int q=0; q<num_lv; q++){
      if(family>1){
         D(q,q,j) = -newlam2(q,j);
      }else{
        D(q,q,j) = 2*newlam2(q,j);
      }
          
    }
  }
  
  //trace of quadratic effect
  for (int i=0; i<n; i++) {
    for (int j=0; j<p;j++){
      if(family>1{
        eta(i,j) -= 0.5*((D.col(j).matrix()*A.col(i).matrix()).trace());
      }else{
      eta(i,j) += ((D.col(j).matrix()*A.col(i).matrix()).trace());}
      
    }
  }
  if(family==0){
    //likelihood
    matrix <Type> e_eta(n,p);
    matrix <Type> B(num_lv,num_lv);
    matrix <Type> v(num_lv,1);
    for (int i=0; i<n; i++) {
      matrix <Type> Q = atomic::matinv(A.col(i).matrix());
      for (int j=0; j<p;j++){
        B = (D.col(j).matrix()+Q);
        v = (newlam.col(j)+Q*u.row(i).transpose());
        Type detB = pow(B.determinant(),-0.5);
        Type detA = pow(A.col(i).matrix().determinant(),-0.5);
        e_eta(i,j) = exp(C(i,j) + 0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
        nll -= y(i,j)*eta(i,j) - e_eta(i,j) - lfactorial(y(i,j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    matrix <Type> Q =  atomic::matinv(A.col(0).matrix());
    //Get lower cholesky.
    matrix <Type> L = Q.llt().matrixL();
    //Transform std normal vector to vector with mean and variance from VA using inverse of cholesky. Not very computationally efficient..
    matrix <Type> znew = u.row(0)*L.inverse();
    REPORT(L);
    REPORT(Q);
    REPORT(znew);
    REPORT(e_eta);
    //input << A.col(0).col(0), A.col(0).col(1);
  }else if(family==1){
    
    // matrix <Type> zetanew(n,p);
    // matrix <Type> B(num_lv,num_lv);
    // matrix <Type> v(num_lv,1);
    matrix <Type> ans(n,p);
    // int length = 2 + num_lv*3 + num_lv*num_lv; //1 for num_lv, 1 for C, q for lambda, lamvda2, u, and q*q for A
    for (int i=0; i<n; i++) {
      //matrix <Type> Q = atomic::matinv(A.col(i).matrix());
      for (int j=0; j<p;j++){
        // vector<Type> inputIntegrand(length);
        // inputIntegrand(0) = num_lv;
        // inputIntegrand(1) = C(0,0);
        // for (int q=0; q<num_lv;q++){
        //   inputIntegrand(2+q) = newlam.col(j)(q);
        //   inputIntegrand(2+num_lv+q) = newlam2.col(j)(q);
        //   inputIntegrand(2+num_lv*2+q) = u.row(i)(q);
        //   for (int q2=0; q2<num_lv;q2++){
        //     matrix <Type> F =A.col(i).matrix();
        //     inputIntegrand(2+num_lv*3+num_lv*q+q2) = F(q,q2);
        //   }
        //   
        // }
        // matrix <Type> newlamtemp(1,num_lv);
        // matrix <Type> newlamtemp2(1,num_lv);
        // matrix <Type> utemp(1,num_lv);
        // for (int q=0; q<num_lv;q++){
        //   newlamtemp(0,q) =  newlam(q,j);
        //   newlamtemp2(0,q) =  newlam2(q,j);
        //   utemp(0,q) = u(i,q);
        // }
        
        // func<Type> f = {newlam.col(j), newlam2.col(j), A.col(i).matrix(),u.row(i),C(i,j),num_lv};
        //ans(i,j) = eval(C(i,j), newlamtemp, newlamtemp2, utemp, A.col(i).matrix());
        //  B = (D.col(j).matrix()+Q);
        //  v = (newlam.col(j)+Q*u.row(i).transpose());
        //  Type detB = pow(B.determinant(),-0.5);
        //  Type detA = pow(A.col(i).matrix().determinant(),-0.5);
        // zetanew(i,j) = iphi(j) + exp(C(i,j) + 0.5*((v.transpose()*atomic::matinv(B)*v).value()-(u.row(i)*Q*u.row(i).transpose()).value()))*detB*detA;
        
        
        //nll -= y(i,j) * eta(i,j) - (y(i,j) + iphi(j))*log(zetanew(i,j)) - iphi(j)*((y(i,j) + iphi(j))/zetanew(i,j)) + lgamma(y(i,j)+iphi(j)) - lfactorial(y(i,j)) + iphi(j)*log(iphi(j)) - lgamma(iphi(j));
        //nll -=  y(i,j)*eta(i,j) - (y(i,j)+iphi(j))*log(iphi(j)+exp(eta(i,j))) + lgamma(y(i,j)+iphi(j)) + iphi(j)*log(iphi(j)) - lgamma(iphi(j)) -lfactorial(y(i,j));
        nll -= -iphi(j)*eta(i,j) - (y(i,j)+iphi(j))*log(1+iphi(j)*ans(i,j)) + lgamma(y(i,j)+iphi(j))+ iphi(j)*log(iphi(j)) - lgamma(iphi(j))  -lfactorial(y(i,j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      REPORT(ans);
    }
  } else if(family==2){
    matrix <Type> mu(n,p);
    for (int i=0; i<n; i++) {
      for (int j=0; j<p;j++){
        mu(i,j) = pnorm(Type(eta(i,j)),Type(0),Type(1));
        nll -= dbinom(y(i,j),Type(1),mu(i,j),true);
        nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    
  } else if(family==3 && zetastruc==1){
    
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
        nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
    
  } else if(family==3 && zetastruc==0){
    
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
        }
        nll -= -0.5*(newlam.col(j)*newlam.col(j).transpose()*A.col(i).matrix()).trace() - (D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*A.col(i).matrix()).trace() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*D.col(j).matrix()*u.row(i).transpose()).value() - 2*(u.row(i)*D.col(j).matrix()*A.col(i).matrix()*newlam.col(j)).value();
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
    }
  }else if(family==4){      
    matrix <Type> ans(n,p);
    for (int i=0; i<n; i++) {
      //operations for integration
      //invert Covariance matrix
      matrix <Type> Q =  atomic::matinv(A.col(i).matrix());
      //Get lower cholesky.
      matrix <Type> L = Q.llt().matrixL();
      //invert it
      matrix <Type> Linv = atomic::matinv(L);
      
      for (int j=0; j<p;j++){
        // func<Type> f = {newlam.col(j), newlam2.col(j), A.col(i).matrix(),u.row(i),C(i,j),num_lv};
        // ans(i,j) = romberg::integrate(f, a, b, n_int, num_lv);//should still catch Infs here..
        // nll -= -lgamma(iphi(j)) + iphi(j)*log(iphi0(j)*y(i,j)) - iphi(j)*eta(i,j) - iphi(j)*y(i,j)* eta2.mean();
        matrix <Type> newlamtemp(1,num_lv);
        matrix <Type> newlamtemp2(1,num_lv);
        matrix <Type> utemp(1,num_lv);
        for (int q=0; q<num_lv;q++){
          newlamtemp(0,q) =  newlam(q,j);
          newlamtemp2(0,q) =  newlam2(q,j);
          utemp(0,q) = u(i,q);
        }
        
        // func<Type> f = {newlam.col(j), newlam2.col(j), A.col(i).matrix(),u.row(i),C(i,j),num_lv};
        ans(i,j) = eval(C(i,j), newlamtemp, newlamtemp2, utemp, A.col(i).matrix(), Linv);//necessary conversion I got from the adaptive_integration example. Otherwise, doesn't compile
        
        nll -=  ( -eta(i,j) - ans(i,j)*y(i,j) )/iphi(j) + log(y(i,j)/iphi(j))/iphi(j) - log(y(i,j)) -lgamma(1/iphi(j));
      }
      // nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      REPORT(ans);
      // Type test = density::MVNORM(A.col(0).matrix())(u.row(0));
      // density::MVNORM_t<Type>nltest(A.col(0).matrix());
      // Type test2 = nltest(u.row(0));
      // REPORT(test);//these were exactly the same..
      // REPORT(test2);
    }
    //report eta2, maybe write the MCI in a function.
  }else if(family==5){
    //here do poisson with numerical integration of exp(eta). Can compare as we have the closed form. Just add a Poisson2 family for now for testing.
    matrix <Type> ans(n,p);
    for (int i=0; i<n; i++) {
      //matrix <Type> Q =  atomic::matinv(A.col(i).matrix());
      //Get lower cholesky.
      // matrix <Type> L = Q.llt().matrixL();
      //invert it
      // matrix <Type> Linv = atomic::matinv(L);
      for (int j=0; j<p;j++){
        // func<Type> f = {newlam.col(j), newlam2.col(j), A.col(i).matrix(),u.row(i),C(i,j),num_lv};
        // ans(i,j) = romberg::integrate(f, a, b, n_int, num_lv);//should still catch Infs here..
        // nll -= -lgamma(iphi(j)) + iphi(j)*log(iphi0(j)*y(i,j)) - iphi(j)*eta(i,j) - iphi(j)*y(i,j)* eta2.mean();
        matrix <Type> newlamtemp(1,num_lv);
        matrix <Type> newlamtemp2(1,num_lv);
        matrix <Type> utemp(1,num_lv);
        for (int q=0; q<num_lv;q++){
          newlamtemp(0,q) =  newlam(q,j);
          newlamtemp2(0,q) =  newlam2(q,j);
          utemp(0,q) = u(i,q);
        }
        
        // func<Type> f = {newlam.col(j), newlam2.col(j), A.col(i).matrix(),u.row(i),C(i,j),num_lv};
        ans(i,j) = eval2(newlamtemp, newlamtemp2, utemp, A.col(i).matrix());//, Linv);//necessary conversion I got from the adaptive_integration example. Otherwise, doesn't compile
        ans(i,j) += C(i,j);
        nll -=  eta(i,j)*y(i,j) - exp(ans(i,j)) -lfactorial(y(i,j));
      }
      nll -= 0.5*(log(Ar(i)) - Ar(i)/pow(sigma,2) - pow(r0(i)/sigma,2))*random(0);
      REPORT(ans);
      // Type test = density::MVNORM(A.col(0).matrix())(u.row(0));
      // density::MVNORM_t<Type>nltest(A.col(0).matrix());
      // Type test2 = nltest(u.row(0));
      // REPORT(test);//these were exactly the same..
      // REPORT(test2);
    }
    
  }
  
  //shrinks LVs, linear ridge
  if(nTol==1){
    for(int q=0; q<(num_lv-1); q++){
      nll += (newlam.row(q).array()*newlam.row(q).array()+newlam2.row(q).array()*newlam2.row(q).array() + newlam.row(q+1).array()*newlam.row(q+1).array()+newlam2.row(q+1).array()*newlam2.row(q+1).array()).sum()*gamma(q);
      nll += (newlam.row(q+1).array()*newlam.row(q+1).array()+newlam2.row(q+1).array()*newlam2.row(q+1).array()).sum()*gamma(q+1);
    }
    
  }else{
    for(int q=0; q<(num_lv-1); q++){
      nll += (newlam.row(q).array()*newlam.row(q).array()+lambda2.row(q).array()*lambda2.row(q).array() + newlam.row(q+1).array()*newlam.row(q+1).array()+lambda2.row(q+1).array()*lambda2.row(q+1).array()).sum()*gamma(q);
      nll += (newlam.row(q+1).array()*newlam.row(q+1).array()+lambda2.row(q+1).array()*lambda2.row(q+1).array()).sum()*gamma(q+1);
    }
  }
  
  //shrinks LVs, quadratic ridge
  for(int j=0; j<lambda2.cols(); j++){
    for(int q=0; q<num_lv; q++){
      nll += lambda2(q,j)*lambda2(q,j)*gamma2(q,j);//should be lamda2...
      //  nll += pow(newlam2(q,j)*newlam2(q,j) + newlam2(q+1,j)*newlam2(q+1,j),0.5)*gamma2(q,j); //should not be hierarchical
      //nll += pow(newlam2(q+1,j)*newlam2(q+1,j),0.5)*gamma2(q+1,j);
    }
  }
  nll -= -0.5*(u.array()*u.array()).sum() - n*log(sigma)*random(0);// -0.5*t(u_i)*u_i
  
  SIMULATE {
    matrix<Type> mu = r0*xr + offset;
    
    if(model<1){
      mu += x*b;
    } else {
      matrix<Type> eta1=x*B;
      int m=0;
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          mu(i,j)+=b(0,j)+eta1(m,0);
          m++;
        }
      }
    }
    
    mu += u*newlam - (u.array()*u.array()).matrix()*newlam2; //intercept(s), linear effect and negative only quadratic term
    matrix<Type>sims(n,p);
    if(family==0){
      for (int j=0; j<p;j++){
        for (int i=0; i<n; i++) {
          sims(i,j) = rpois(exp(mu(i,j)));
        }
      }
    }
    // }else if(family==1){
    //
    // }
    
    REPORT(sims);          // Report the simulation
    
  }
  return nll;
}

//for (int q=0; q<num_lv;q++){//
//  integrand<Type> f = {newlam(q,j), newlam2(q,j)};
//ans(i,j) += log(gauss_kronrod::integrate(f, Type(-5),Type(5)));
