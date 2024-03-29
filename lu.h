/******************************************************************
*  三对角线性方程组的基本求解方法程序
*　　－－－LU求解算法函数定义的头文件
*  
*  文件名：lu.h   
*  版    权：共享
*
*  算法编写：刘智翔   
*  程序编写：郑汉垣     
*  创建时间：2010.12.23 01:00 
*  版本号：Ｖ1.５
*  修改人：郑汉垣
*  修改时间：2010.12.24 0９:４0
*******************************************************************/

/* 避免重复包含该头文件 */
#ifndef _LU_H_
#define _LU_H_

/*包含三对角线性方程组求解算法中使用的宏代换定义的公用头文件*/
#include "trieq.h"

/****************************************************************
*  函数功能：三对角线性方程组的正向、反向、双向ＬＵ分解法求解				*
*  输入参数：A(N,4)是A是由向量a,b,c,d组成的矩阵，N是向量的长度			*
*  输出参数：线性方程组的解Ｘ(N), 函数值为０则正常求解					*
*  注：在求解过程中，A(N,4)已经被改变,是因为节省内存的原因				*
*****************************************************************/
int LinearEqRecduce(double *A,long n, double *Mat2);  //对Matq(q,4)对应的三对角方程进行约化
int forwardLU(double *A,long n, double *X);      //主进程求解2*p行数的缩减方程，获得2*p个解，x2p(2*p)
int bilateralLU(double *A,long n, double *X);      //主进程求解2*p行数的缩减方程，获得2*p个解，x2p(2*p)
int LinearEqBacktrack(double *A, long n, double *X2, double *X);   //对Mat(n,4)对应的三对角方程进行回代
 
#endif   //_LU_H_
