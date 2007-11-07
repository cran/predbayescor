#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <values.h>
/* this function calculates the joint distribution of a feature */
/* here only calculate the part involved with theta */
double urn1(double alpha, double theta,int I0,int O0,int I1,int O1)
{ if(theta==0){
     if(I0!=0 || I1!=0) return -MAXDOUBLE;
     else return 
        lgammafn( alpha  + O0 ) +
        lgammafn( alpha  + O1 ) -
        2 * lgammafn( alpha  );
  }
  
  if(theta==1){
     if(O0!=0 || O1!=0) return -MAXDOUBLE;
     else return
        lgammafn( alpha  + I0 ) +
        lgammafn( alpha  + I1 ) -
        2 * lgammafn( alpha  );
  }   
  return 
    lgammafn( alpha * theta + I1 ) + 
    lgammafn(alpha * theta +I0) +
    lgammafn( alpha * (1-theta) + O0 ) +
    lgammafn( alpha * (1-theta) + O1 ) -
    2 * lgammafn( alpha * theta ) - 
    2 * lgammafn( alpha * ( 1-theta ) ) ; 
}
/*this function gives the log probability of Bernoulli distribution */
double bern(int x, double phi)
{  if(x==1) return log(phi); 
   else return log(1-phi); 
}
/*this function finds the maximum value of a vector */
double findmax(int la,double a[])
{ int i, id_max=0;
  for(i=1;i<la;i++)
  {  if(a[i]>a[id_max]) id_max = i;
  }
  return a[id_max];
}
/*this function calculates the log sums of a vector expressed in log form */
double log_sum_exp(int la, double a[])
{ double m,s=0;
  int i;
  m = findmax(la,a);
  for(i=0;i<la;i++) s += exp(a[i]-m);
  return log(s) + m;
}
/*this function transforms a value "x" in (0,1) to another
value in (0,1) such that the new one is closer to "theta0". 
"n" is used to control how close, the bigger the closer */
double tfm_unif(double x,double theta0,int n)
{  double a,aPN;
   int N;
   N = 2*n+1;
   a = 1.0 / (1.0 + R_pow(1.0/theta0 - 1, 1.0/N) );
   aPN=R_pow_di(a,N);
   return (R_pow_di(x-a,N)+ aPN)/(R_pow_di(1-a,N)+ aPN);
}

/*the inverse transformation of tfm_unif */
double tfm_unif_inv(double theta,double theta0,int n)
{  double a,aPN,tmp,abs_tmp;
   int N;
   N = 2*n+1;
   a = 1.0 / (1.0 + R_pow(1.0/theta0-1, 1.0/N) );
   aPN = R_pow_di(a,N);
   tmp = theta*(R_pow_di(1-a,N)) + aPN*(theta - 1);
   if(tmp>0)
       return R_pow(tmp,1.0/N) + a;
   else
       return -R_pow(-tmp,1.0/N) + a;
}

/*the derivative of "tfm_unif" with respect to x */
double tfm_unif_dir (double x,double theta0, int n)
{  double a;
   int N;
   N = 2*n+1;
   a = 1.0 / (1.0 + R_pow(1.0/theta0 - 1, 1.0/N) );
   return (R_pow_di(x-a,N-1)*N)/(R_pow_di(1-a,N)+ R_pow_di(a,N));
}
/******************************************************************************/
/*this function does the prediction*/
/*the notations of arguments are almost consistent with the notations in paper*/
/*
k : number of features retained
no_test: number of test cases
no_train0: number of training cases of class 0
no_train1: number of training cases of class 1
no_alpha:  number of alpha's used 
no_theta: number of theta's used in numerical integration
test: test dataset 
theta_range: the starting and ending point of the prior of theta
I1: the vector of the number of 1's for each feature in training data of class 1
I0: the vector of the number of 0's for each feature in training data of class 0
log_prior: the vector of log prior probabilities over the alpha's, which may
           have been adjusted by information from selection
alpha_bank: the vector of all alpha's 
prob_pred: the probabilities of class 1 for each test case, of length no_test[0]
no_theta2: 2*no_theta[0]+1, for convenient usage in simpson rule integration
log_post: the 3-dimensional array of the log joint distribution of training data
          at each alpha (1st dim), each theta (2nd dim),and each feature (3rd dim)
theta: the 2-dimensional array of all theta's of all features (2nd dim)
log_tfm_dir: the 3-dimensional array of the log derivative of the transformation
phi0 : the Bayesian predictive means, a 3-dimensional array at each alpha (1st dim), 
       each theta (2nd dim),and each feature (3rd dim), for class 0
phi1: similar as phi0, but for class 1
dg_tfm: the power of transformation, could be non-negative integer
post-show: the indicator whether to calculate the marginal posterior
           distribution of alpha's      
*/
 
void predict_int(  int* k, int* no_test, int* no_train0, int* no_train1,
                  int* no_alpha,int* no_theta, int test[no_test[0]][k[0]], 
		  double theta_range[2],int I1[k[0]],int I0[k[0]], 
		  double log_prior[no_alpha[0]], double alpha_bank[no_alpha[0]],
		  double prob_pred[no_test[0]], int* no_theta2,
		  double log_post[no_alpha[0]][no_theta2[0]][k[0]],
		  double theta[][k[0]],double log_tfm_dir[][k[0]],
		  double phi0[no_alpha[0]][no_theta2[0]][k[0]],
		  double phi1[no_alpha[0]][no_theta2[0]][k[0]],
		  int* dg_tfm,double *log_post_alpha)
{   int a,i,j,b, no_train=no_train0[0]+no_train1[0];
    double h,xtheta_b,xtheta_e,*log_prob_theta0,*log_prob_theta1,
          *log_prob_alpha0, *log_prob_alpha1,*urn_p2,logp0,logp1,mean_theta,
          logtwo = log(2),log_sum_post;
    /*allocate memory for log_prob_theta0 --- a vector of log predictive
       probabilites for a feature over different theta's, which is reused for
       all features. It is also used to hold the log posterior probabilities
       of a feature over different theta's */
    log_prob_theta0 = (double *)R_alloc(no_theta2[0],sizeof(double));
    log_prob_theta1 = (double *)R_alloc(no_theta2[0],sizeof(double));
    /*allocate memory for log_prob_alpha0 --- vector of log predictive
       probabilites for a test case over different alpha's, which is reused 
       for all test cases. */
    log_prob_alpha0 = (double *)R_alloc(no_alpha[0],sizeof(double));
    log_prob_alpha1 = (double *)R_alloc(no_alpha[0],sizeof(double));
    urn_p2 = (double *)R_alloc(no_alpha[0],sizeof(double));
    
    /*listing all thetas used in integration by transforming from the evenly
       spaced points between xtheta_b --- the inverse transformation of 
       theta_range[0]
       and xtheta_e --- the inverse transformation of theta_range[1]. The
       transformation makes the theta's closer to the mean of feature j */
    for(j=0;j< k[0];j++)
    {  /*the mean of feature j*/
       mean_theta = (double)(I0[j]+I1[j])/no_train;
       xtheta_b = tfm_unif_inv(theta_range[0],mean_theta,dg_tfm[0]);
       xtheta_e = tfm_unif_inv(theta_range[1],mean_theta,dg_tfm[0]);			
       h = (double)(xtheta_e- xtheta_b)/(2*no_theta[0]); 
       
       for(i=0;i<*no_theta2;i++)
       {  /*tranformating values btw xtheta_b, xtheta_e for 
               making up thetas[i][j]  */
          theta[i][j] = tfm_unif(xtheta_b+i*h,mean_theta,dg_tfm[0]);
          /*calculating the log derivative at thetas[i][j] */ 
	  log_tfm_dir[i][j] = log(tfm_unif_dir(xtheta_b+i*h,mean_theta,
	                          dg_tfm[0]));	  
       }
    }
    /*calculate the second part of log joint (post) distribution which 
      is without theta */
    for(a=0;a< no_alpha[0]; a++)
    {  urn_p2[a]= 2 * lgammafn(alpha_bank[a]) - 
                 lgammafn(no_train0[0] + alpha_bank[a]) - 
                 lgammafn(no_train1[0] + alpha_bank[a]);
    }
    /* calculate the phi's and post of phi on the grid of thetas and caculate 
       the posterior of alpha using simpson rule if it is requried*/
    for(a=0;a< no_alpha[0]; a++)			      
    {  log_post_alpha[a] = log_prior[a] + k[0]*urn_p2[a];
       for(j=0;j< k[0];j++)
       {  for(i=0;i<*no_theta2;i++)
          {   /*the predictive means of feature j at alpha.bank[a] 
                  and theta[i][j], respectively for class 0 and class 1 */
              phi0[a][i][j] = (I0[j] + alpha_bank[a]*theta[i][j]) /
	                      (alpha_bank[a] + no_train0[0]);
	      phi1[a][i][j] = (I1[j] + alpha_bank[a]*theta[i][j]) /
	                      (alpha_bank[a] + no_train1[0]);
	      /*the post disctribution at alpha[a], theta[i][j]*/
	      log_post[a][i][j] = urn1(alpha_bank[a],theta[i][j],
				   I0[j],no_train0[0]-I0[j],
				   I1[j],no_train1[0]-I1[j]) + 
			           log_tfm_dir[i][j];
	       
	      log_prob_theta0[i] = log_post[a][i][j];
	      if(i!= 0 & i!= no_theta2[0]-1)
	      {  log_prob_theta0[i] += logtwo * (i-i/2*2+1);
	      }
	  }       
          log_post_alpha[a] += log_sum_exp(no_theta2[0],log_prob_theta0);
       }
    }
    /*gives the log sum of posterior distribution of alpha */
    log_sum_post = log_sum_exp(no_alpha[0],log_post_alpha);
    /*standardize the posterior of alpha by substracting the log sum */
    for(a=0;a< no_alpha[0]; a++)
       log_post_alpha[a] -= log_sum_post;
        
    /*start making prediction */ 
    /*iterated for test cases */
    for(b=0;b < no_test[0];b++)
    {   /*iterated for alpha's */
        for(a=0;a < no_alpha[0]; a++)
        {  log_prob_alpha0[a]=log_prior[a] + k[0]*urn_p2[a];
	   log_prob_alpha1[a] = log_prob_alpha0[a];
	   /*iterated for features */
	   for(j=0;j< k[0]; j++)
	   {  /* calculate the log predictive probabilities for feature j at
	         theta[i][j], for later use in the integration with 
	         simpson rule */   
	      for(i=0;i< no_theta2[0];i++)
	      {   log_prob_theta0[i] = 
	              bern(test[b][j],phi0[a][i][j])+log_post[a][i][j];
		  log_prob_theta1[i] = 
		      bern(test[b][j],phi1[a][i][j])+log_post[a][i][j];
		  if(i!= 0 & i!= no_theta2[0]-1)
		  {   log_prob_theta0[i] += logtwo * (i-i/2*2+1);
		      log_prob_theta1[i] += logtwo * (i-i/2*2+1);
		      
		  }
	      }       
              /*calculate the log sum of predictive probabilities for 
                feature j. i.e. doing simpson rule integration, and  
                added to log predicitive prob at alpha[a]*/
	      log_prob_alpha0[a]+=log_sum_exp(no_theta2[0],log_prob_theta0);     
	      log_prob_alpha1[a]+=log_sum_exp(no_theta2[0],log_prob_theta1);
	   }	   	   
	 }
	 /*calculate the log sum of predictive probabilites acrossing alpha */
	 logp0=log_sum_exp(no_alpha[0],log_prob_alpha0)+log(no_train0[0]+1);
	 logp1=log_sum_exp(no_alpha[0],log_prob_alpha1)+log(no_train1[0]+1);
	 /*calculate the predictive probability for case b */
	 prob_pred[b] = 1.0/(1+exp(logp0-logp1)); 
    }    
}
/*this function is to be called by R. It only allocates memory for
   phi0,phi1,theta, log_post, log_tfm_dir which will be then passed to
   "predict_int". This makes it convenient to refer to the above
   array element after defining dimensions by a new function. 
   Please go to predict_int 
   for explanation of arguments */		    
void predict_int_c( int* k, int* no_test, int* no_train0, int* no_train1,
                   int* no_alpha,int* no_theta, int test[no_test[0]][k[0]],
		   double theta_range[2],int I1[k[0]], 
		   int I0[k[0]], double log_prior[no_alpha[0]],
		   double alpha_bank[no_alpha[0]],double prob_pred[no_test[0]],
		   int *dg_tfm,double *log_post_alpha)
{   double *log_post,*theta,*log_tfm_dir,*phi0,*phi1;
    int no_theta2=2*(no_theta[0])+1;
    log_post = (double *) R_alloc(no_alpha[0]*no_theta2*k[0],sizeof(double));
    phi0 = (double *) R_alloc(no_alpha[0]*no_theta2*k[0],sizeof(double));
    phi1 = (double *) R_alloc(no_alpha[0]*no_theta2*k[0],sizeof(double));
    theta = (double *) R_alloc(no_theta2*k[0],sizeof(double));
    log_tfm_dir = (double *) R_alloc(no_theta2*k[0],sizeof(double));
    predict_int(k,no_test,no_train0,no_train1,no_alpha,no_theta,test,
                theta_range,I1, I0,log_prior,alpha_bank,prob_pred,&no_theta2,
		log_post,theta,log_tfm_dir,phi0,phi1,dg_tfm,log_post_alpha);      
}
/****************************************************************************/
/*this fucntion calculates the probability of a feature has low correlation
 given alpha and theta */ 
double prob_lowcor_theta(int n0[1], int n1[1],int no_I0[1], int I0[],int b[1], 
                         int I0_max[1],double alpha[1],double theta[1])
{  double prob1[no_I0[0]],prob0[I0_max[0]+1],cprob0[I0_max[0]+1],a0,a1,
          prob_highcor,lprob_alphatheta;
   int i;
   if(theta[0]==0 || theta[0]==1) return 1;
   a1 = alpha[0]*theta[0];
   a0= alpha[0]-a1;
   lprob_alphatheta = lgammafn(alpha[0])-lgammafn(a1)-lgammafn(a0);
   prob1[0] = exp(lchoose(n1[0],b[0])+lgammafn(a1+b[0])+
                  lgammafn(a0+n1[0]-b[0])-lgammafn(alpha[0]+n1[0])+
	          lprob_alphatheta);
   prob0[0] = exp(lgammafn(a1)+lgammafn(a0+n0[0])-
                  lgammafn(alpha[0]+n0[0]) + lprob_alphatheta);
   cprob0[0] = prob0[0];   
   for(i=0; i < I0_max[0]; i++)
   {  prob0[i+1]=(double) prob0[i]*(a1+i)/(i+1)*(n0[0]-i)/(a0+n0[0]-i-1);
      cprob0[i+1] = cprob0[i] + prob0[i+1];
      /*Rprintf("%f\n",cprob0[i]);*/
   }
   prob_highcor = prob1[0] * cprob0[I0[0]];
   for(i=0; i < no_I0[0]-1; i++)
   {  prob1[i+1]=(double) prob1[i]*(a1+b[0]+i)/(b[0]+i+1)*(n1[0]-(b[0]+i))/
                 (a0+n1[0]-(b[0]+i)-1);
      prob_highcor += prob1[i+1] * cprob0[I0[i+1]];
      /* Rprintf("%f\n",prob1[i]);*/
   }
   /* Rprintf("%f,%f,%f\n",alpha[0],theta[0],prob_highcor);*/
   return 1-2*prob_highcor;        
}
/*this function calculates that a feature has low correlation with response 
given alpha */		
double prob_lowcor_alpha(int no_theta[1], int n0[1], int n1[1], 
                        int no_alpha[1], int no_I0[1], int I0[],
		        int b[1], int I0_max[1],double alpha[1], 
		        double theta_range[2])
{   double prob,h,thetas[2*no_theta[0]+1],thetas_dir[2*no_theta[0]+1],
           xtheta_b,xtheta_e;
    int i;
    xtheta_b = theta_range[0]; /*0.5-R_pow((0.5-theta_range[0])/4,1.0/3);*/
    xtheta_e = theta_range[1]; /*0.5+R_pow((theta_range[1]-0.5)/4,1.0/3);*/
    h = (double) (xtheta_e - xtheta_b)/(2*no_theta[0]);
    for(i=0;i<2*no_theta[0]+1;i++)
    {  thetas[i] = xtheta_b+i*h; /*R_pow_di(xtheta_b+i*h-0.5,3)*4 + 0.5;*/
       thetas_dir[i] = 1; /*12 * R_pow_di(xtheta_b+i*h-0.5,2);*/
    }
    /* doing simpson rule integration */
    prob = prob_lowcor_theta(n0,n1,no_I0,I0,b,I0_max,alpha,&thetas[0]) 
           * thetas_dir[0];
    prob += prob_lowcor_theta(n0,n1,no_I0,I0,b,I0_max,alpha,
	   &thetas[2*no_theta[0]]) * thetas_dir[2*no_theta[0]];
    	   
    for(i=1;i< no_theta[0]+1;i++)
    {   prob += 4*prob_lowcor_theta(n0,n1,no_I0,I0,b,I0_max,alpha,
                   &thetas[2*i-1]) * thetas_dir[2*i-1];
    }
    for(i=1;i< no_theta[0];i++)
    {	prob += 2*prob_lowcor_theta(n0,n1,no_I0,I0,b,I0_max,alpha,
                   &thetas[2*i]) * thetas_dir[2*i];
    }
    return prob*h/3/(theta_range[1]-theta_range[0]);    
}
/*this function is to be called by R*/		
void prob_lowcor(int no_theta[1], int n0[1], int n1[1], 
                 int no_alpha[1], int no_I0[1], int I0[],
		 int b[1],int I0_max[1], double alpha_bank[], 
		 double theta_range[2], double prob_lowcor[no_alpha[0]])
{  int i;
   for(i=0;i< no_alpha[0];i++)
   {  prob_lowcor[i] = 
         prob_lowcor_alpha(no_theta,n0,n1,no_alpha, 
	                   no_I0,I0,b,I0_max,&alpha_bank[i],theta_range);
   }
}		  
/*****************************************************************************/
/* this function calculates the correlation */
double calcorl(int I0,int I1,int n, double y_bar, double y_std)
{  if(I0+I1 == 0 || I0+I1 == n) return 0;
   else
      return (-y_bar*I0+(1-y_bar)*I1) / 
             ( y_std * R_pow(I0+I1-R_pow_di(I0+I1,2)/n,0.5)) ;
}
/* this function determines the bounds of I0 which makes high correlation */
void select_bd(int n0[1],int n1[1],double gamma_sel[1],int b[1], int I0[])
{  int i,n;
   double y_bar,y_std;
   n = n0[0] + n1[0];
   y_bar = (double) n1[0] / n;
   y_std = R_pow(n1[0]*(1-y_bar),0.5); 
   I0[0]=0;   
   for(;;)
   {  if(calcorl(I0[0]+1,b[0],n,y_bar,y_std) < gamma_sel[0]) break;
      else I0[0]++;
   }
   for(i=1; i < n1[0]-b[0]+1; i++)
   {  I0[i] = I0[i-1];
      for(;;)
      {  if(calcorl(I0[i]+1,b[0]+i,n,y_bar,y_std) < gamma_sel[0]) break;
         else I0[i]++;
      }  
   }
}
/*****************************************************************************/
/* only for debugging in R */
void tfm_unif2 (double *x,double *theta0, int *n)
{  Rprintf("%f\n",tfm_unif(x[0],theta0[0],n[0]));}
/* only for debugging in R */
void tfm_unif_inv2 (double *theta,double *theta0, int *n)
{  Rprintf("%f\n",tfm_unif_inv(theta[0],theta0[0],n[0]));}
/* only for debugging in R */
void tfm_unif_dir2 (double *x,double *theta0, int *n)
{  Rprintf("%f\n",tfm_unif_dir(x[0],theta0[0],n[0]));}
