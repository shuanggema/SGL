#include "stdio.h"
#include "math.h"
#include "stdlib.h"

int Sum(int *G, int n)
{
  int i,value=0;
  for( i = 0; i < n; i++)
     value += G[i];
     
  return value;
}

double Norm(double *z,int n)
{
  int i;
  double sum=0,value;
  for ( i = 0; i < n; i++)
  {
    sum += z[i]*z[i];
  }
  value = pow(sum,0.5);
  return value;
}

double Norm_b(double *Beta, int *G, int pos_s, int j)
{
  double sum = 0, value;
  int jj;
  for ( jj = 0; jj< G[j]; jj++ )
  {
     sum += Beta[pos_s+jj]*Beta[pos_s+jj];
  }
  value = pow(sum,0.5);
  return value;
}

double Inner_prod(double *X,double *Beta,int totalp)
{
  int j;
  double sum=0;
  
  for (j = 0; j < totalp; j++)
  {
    sum += X[j]*Beta[j];
  }
  return sum;
}

void S_func(double *z,double *X,double *Beta,double *r,double a, double c,
            int pos_s,int pos_e,int d_j,int n,int totalp)
{
  int i,j;
  double norm, cons, diff[d_j],sum;
 
  norm = Norm(z,d_j);
  cons = (1 - c/norm);
  
  if (cons < 0)
    cons=0;
  
  for ( j = pos_s; j <= pos_e; j++)
  {
    diff[j-pos_s] = cons*z[j-pos_s]/a - Beta[j]; 
    Beta[j] = cons*z[j-pos_s]/a;
  }
  //update r here;
  for ( i = 0; i < n; i++)
  {
    sum = 0;
    for ( j = pos_s; j <= pos_e; j++)
    {
      sum += X[totalp*i+j]*diff[j-pos_s];
    }
    r[i] -= sum;
  }
}

void Update_beta(double *X, double *r, double *W2, double *Beta, double lambda1, double lambda2,
                 int n,int totalp,int *pos_s,int *pos_e,int *G, int J, int P, int D)
{
  int i,j,s=pos_s[J],e=pos_e[J];
  double z[G[J]],a,c;
  
  for ( j = s;j <= e;j++)
  {
    //let i-pos_s be the position for z[j];
    z[j-s] = 0;
    for ( i = 0; i < n; i++)
    {
      z[j-s] += X[totalp*i+j]*r[i]/n;
    }
    z[j-s] += Beta[j];
  }
  if ( J == 0 )
  {
     a = 1 + lambda2*D*W2[J]/G[J];
     c = pow(G[J],0.5)*lambda1-lambda2*D/pow(G[J],0.5)*W2[J]*Norm_b(Beta,G,pos_s[J+1],J+1)/pow(G[J+1],0.5);
  }    
  else if ( J == P-1 )
  {
     a = 1 + lambda2*D*W2[J-1]/G[J];
     c = pow(G[J],0.5)*lambda1-lambda2*D/pow(G[J],0.5)*W2[J-1]*Norm_b(Beta,G,pos_s[J-1],J-1)/pow(G[J-1],0.5);
  }
  else 
  {
     a = 1 + lambda2*D*(W2[J]+W2[J-1])/G[J];
     c = pow(G[J],0.5)*lambda1-lambda2*D/pow(G[J],0.5)*(W2[J-1]*Norm_b(Beta,G,pos_s[J-1],J-1)/pow(G[J-1],0.5)+W2[J]*Norm_b(Beta,G,pos_s[J+1],J+1)/pow(G[J+1],0.5));
  }
    
    
  S_func(z,X,Beta,r,a,c,s,e,G[J],n,totalp);
}

double Gslasso_sp(double *Y, double *X, int *G, int *Param, double *Lambda1,double *Lambda2, double *W2, double *Epsilon,
              double *Beta)
{
  int i,j,k,pos_s[Param[1]],pos_e[Param[1]],count,n=Param[0],p=Param[1],totalp=Param[2],iter=Param[3],nlambda=Param[4],d=Param[5];
  //p is the number of group, totalp is the number of covariates;
  double **Xt,r[Param[0]],old_beta[Param[2]], diff, diff_beta[Param[2]],**beta,lambda1,lambda2;

  Xt = malloc(sizeof(double *) * Param[0]); 
  for ( i = 0; i < Param[0]; i++)   
  { 
     Xt[i] = malloc(sizeof(double) * Param[2]); 
  }
    
  for ( i = 0; i < Param[0]; i++ )   
  {
    for ( j = 0; j < Param[2]; j++ )   
    {
       Xt[i][j] = *(X+Param[2]*i+j); 
    }
  }   
  
  beta = malloc(sizeof(double *) * nlambda); 
  for ( k = 0; k < nlambda; k++)   
  { 
     beta[k] = malloc(sizeof(double) * totalp); 
  }

  for ( j = 0; j < totalp; j++ )   
  {
    beta[0][j] = 0; 
  }
  
  for( j = 0; j < p; j++)
  {
    if ( j == 0)
       pos_s[j]=0;
    else
       pos_s[j]=Sum(G,j);
    pos_e[j]=pos_s[j]+G[j]-1;
  }

  for ( k = 0; k < nlambda; k++)
  {
    if ( k != 0 )
    {
      for ( j = 0; j < totalp; j++ )
      {
        beta[k][j] = beta[k-1][j];
      }
    } 

    // do initialization in outer loop;
    for( i = 0; i < n; i++)
    {
      r[i] = Y[i] - Inner_prod(Xt[i],beta[k],totalp);
    }
  
    for ( j = 0; j < p; j++)
    {
      old_beta[j] = beta[k][j];
    }     
    
    lambda1 = Lambda1[k];
    lambda2 = Lambda2[k];    
    count = 0;
    do
    { 
    
      for ( j = 0; j < p; j++)
      {
        Update_beta(X, r, W2, beta[k], lambda1,lambda2, n, totalp, pos_s, pos_e, G,j,p,d);     
      }

      for ( j = 0; j < totalp; j++)
        diff_beta[j]=old_beta[j]-beta[k][j];
      
      diff = Norm(diff_beta,totalp);        
      count=count+1;

      for ( j = 0; j < totalp; j++)
      {
        old_beta[j] = beta[k][j];
      }       
      
    }
    while ( count <= iter && diff > Epsilon[0]); 
    
    printf("iter=%d\n",count); 
    
  }
  
  for ( k = 0; k < nlambda ; k++ )
    for ( j = 0; j < totalp; j++ )
      Beta[k*totalp+j]= beta[k][j];

  for ( i = 0 ; i < n; i++)
  {  
    free(Xt[i]);
  }  
  free(Xt);
  for ( k = 0 ; k < nlambda; k++) 
  {
    free(beta[k]);
  }
  free(beta);
  return 0;
}


