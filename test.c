/******************************************************************
*  三对角线性方程组的基本求解方法程序
*　　－－－主函数中的测试函数定义文件
*  
*  文件名：test.c   
*  版    权：共享
*
*  程序编写：郑汉垣     
*  创建时间：2010.12.23 01:00 
*  版本号：Ｖ1.3
*  修改人：郑汉垣
*  修改时间：2010.05.17 08:30
*******************************************************************/

#include "test.h"

/****************************************
*   以下是主程序中测试中的子函数	*　　  
 * *************************************/
/* 准备数据A与Ｎ和Ｍ*/

int preDate(double *A,long n,int myid, int numprocs)
{
/* 开始准备数据..... */
   long i,l;
   //  A=(double *)malloc(doublesize*n*4);
    for(i = 1; i <= n; i ++)
    {
      a(i)=-1.0;
      b(i)= 2.0;
      c(i)=-1.0;
      d(i)= 0.0;
    }
   if(myid==0)
   { 
      a(1)= 0.0;
      d(1)= 1.0;
   }
   if(myid==numprocs-1)
   {
      c(n)= 0.0;
      d(n)= 1.0;
   }
    return 0;
}

/*生成新的矩阵Ｂ，并将源矩阵Ａ的值复制到新生成的矩阵Ｂ中　*/
double *copyMat(double *A, long n, long m)
{   
	long i,j;
	double *B;
	
   B=(double *)malloc(doublesize*n*4);	
    
    for(i = 0; i < n; i ++)
    {
      for(j = 0; j < m; j ++)
      {
          *(B+i*m+j) = *(A+i*m+j);
       }
    }
    
    return B;   /* 向量的数据个数 */
}

/* 打印，测试读入的数据是否正确 */    
 int  printMat(double *A,long n,int myid, int numprocs)
 {   
    long i;

    printf("=======myid=%d=======\n",myid); 
    printf("\nMatrix A:\n a \t\t b \t\t c \t\t d \n");

     /* 打印a(i),b(i),c(i),d(i) */    
    for(i=1;i<=2;i++)
    {
          printf("%lf\t %lf\t %lf\t %lf \n", a(i),b(i),c(i),d(i));
    }
    printf("..............................\n");
    i=n;
    printf("%lf\t %lf\t %lf\t %lf \n", a(i),b(i),c(i),d(i));
    printf("\n");
   
    return 0;
}

/* 输出求出的三对角线性方程组的解 */
int  printVal(double *X,long n,int myid, int numprocs)
{
    long i;

    printf("=======myid=%d=======\n",myid); 
    printf("X=[");
    for(i=1;i<=2;i++)
    {
      printf("%lf,",x(i));
    }   
      
    printf("......,%lf]\n\n",x(n));
    return 0;
}

/*  回收分配的内存 */
int  recovery(double *A,double *X)
{
    free(A);
    free(X);
    return 0;
}

/* 程序文本结束 */
