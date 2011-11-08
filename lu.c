﻿/******************************************************************
*  三对角线性方程组的基本求解方法程序
*　　－－－LU求解算法函数定义文件
*  
*  文件名：lu.c   
*  版    权：共享
*
*  算法编写：刘智翔   
*  程序编写：郑汉垣     
*  创建时间：2010.12.23 01:00 
*  版本号：Ｖ2.0
*  修改人：郑汉垣
*  修改时间：2011.05.18 01:30
*******************************************************************/

#include "lu.h"

/*===============以下是主进程和从进程调用的三个过程===============================
  LinearEqRecduce(A, N, Mat2) 对子块方程进行约化，获取两行的子缩减方程
  forwardLU(A,N,X)            三对角线性方程组正向LU分解法求解
  LinearEqBacktrack(A, N, X2, X) 
             给定子块方程的解x2(x1,xn)，对子块方程进行回代，获取子块方程的解x(n)
================================================================================= */

/*****************************************************************
  函数功能：对子块方程进行约化，获取两行的子缩减方程(8个数据)
  输入参数：A是由向量a,b,c组成的矩阵，右向量为d, n是向量的长度
  输出参数：线性方程组的解Mat2(n,4), 函数值为０则正常求解
  注：在求解过程中，a,b,c,d,Mat(n,4)都已被改变
******************************************************************/
int LinearEqRecduce(double *A,long n, double *Mat2)  //对Matq(q,4)对应的三对角方程进行约化
{

   double  t;
   long  i, j; 
   
   //自上而下的消元（追）
   for(i=3; i<=n; i++)
   {
         t = -a(i)/b(i-1);
      a(i) = a(i-1) * t;
      b(i) = b(i) + c(i-1)*t;
      d(i) = d(i) + d(i-1)*t;
   }

   //自下而上的消元（赶）
   for(i=n-2; i>=2; i--)
   {
         t = -c(i)/b(i+1);
      c(i) = c(i+1) * t;
      a(i) = a(i) + a(i+1)*t;
      d(i) = d(i) + d(i+1)*t;
   }
         i = 1 ;
         t = -c(i)/b(i+1);
      c(i) = c(i+1) * t;
      b(i) = b(i) + a(i+1)*t;
      d(i) = d(i) + d(i+1)*t;

  //返回Mat2  子块方程的缩减方程：[a1,b1,c1,d1;an,bn,cn,dn]
    Mat2[0] = a(1);
    Mat2[1] = b(1);
    Mat2[2] = c(1);
    Mat2[3] = d(1);

    Mat2[4] = a(n);
    Mat2[5] = b(n);
    Mat2[6] = c(n);
    Mat2[7] = d(n);

    return 0;
}

/****************************************************************
  函数功能：三对角线性方程组正向LU分解法求解
  输入参数：A是由向量a,b,c组成的矩阵，右向量为d, n是向量的长度
  输出参数：线性方程组的解x(n), 函数值为０则正常求解
  注：在求解过程中，a,b,c,d,x都已被改变
****************************************************************/
int forwardLU(double *A,long n, double *X)    //主进程求解2*p行数的缩减方程，获得2*p个解，x2p(2*p)
{
    long   i;       //中间变量 
    double t;
         
    //自上而下的消元
    c(1) = c(1) / b(1);
    d(1) = d(1) / b(1);
    for( i = 2; i<=n-1; i++)
    {
        t = b(i)-(a(i)*c(i-1));
        c(i) = c(i)/t;
        d(i) = (d(i)-(a(i)*d(i-1)))/t;
     }
        i = n;
        t = b(i)-(a(i)*c(i-1));
        d(i) = (d(i)-(a(i)*d(i-1)))/t;

    //自下而上的回代    
    x(n) = d(n);
    for(i=n-1; i>=1; i--)
    {
       x(i) = d(i)-c(i)*x(i+1);
    }

   //返回三对角线性方程组的解x(n),但是a,b,c,d向量数据都被更改了。
   return 0; 
}

/********************************************************************
  函数功能：给定子块方程的解x2(x1,xn)，对子块方程进行回代，获取子块方程的解x(n)
  输入参数：A是由向量a,b,c组成的矩阵，右向量为d, n是向量的长度,x2是子块方程的解
  输出参数：线性方程组的解x(n), 函数值为０则正常求解		
  注：在求解过程中，a,b,c,d,x都被改变	         	
********************************************************************/
int LinearEqBacktrack(double *A, long n, double *X2, double *X)   //对Mat(n,4)对应的三对角方程进行回代
{
   long i;

   x(1)=X2[0];
   x(n)=X2[1];

   for(i=2;i<=n-1;i++)
   {
      x(i)= ( d(i)-a(i)*x(1)-c(i)*x(n) )/b(i);
   }
  return 0;
}
