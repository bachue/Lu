/******************************************************************
*  三对角线性方程组的基本求解方法程序
*　　－－求解算法函数的主函数定义文件
*  
*  文件名：main.c   
*  版    权：共享
*
*  程序编写：郑汉垣     
*  创建时间：2010.12.23 01:00 
*  版本号：Ｖ1.5
*  修改人：郑汉垣
*  修改时间：2011.05.17 22:40
*******************************************************************/

/* 包含主函数中的测试函数定义的头文件 */
#include "test.h"

/* 包含LU求解算法函数定义的头文件 */
#include "lu.h"

#include <stdio.h>
#include <math.h>

//  #include <omp.h>

#include "mpi.h"

int main(int argc,char *argv[])
{
    int myid,numprocs;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int nthreads, tid;
    char buf[32];

    long m,p,q,k,r;   //线性方程总个数，进程数，子块行数，线程数，子子块行数
    double  starttime,endtime,walltime;
    
    double *A;   /* 按行方式存储的三对角矩阵与右向量(计算用) */
    double *X;   /* 求值向量 */
    
    double *Mat2p;   /* 按行方式存储的三对角矩阵与右向量(计算用) */
    double *X2p;   /* 求值向量 */
     
    double *Mat2;   /* 按行方式存储的三对角矩阵与右向量(计算用) */
    double *X2;   /* 求值向量 */
    double *Mat;

    long i,j,n;
    int l,nn,ret,root,count;

 
   //m=1024*1024*64;   //m=2**n  (n=26)
   nn=28;
   m=1;
   for(k=1;k<=nn;k++) 
   {
	m=2*m;
    }
   printf("m=%ld\n",m);

    MPI_Init(&argc,&argv);   /*MPI init */
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);     /* processor id */
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* processor size */
    MPI_Get_processor_name(processor_name,&namelen); 

    p=numprocs;  //进程总数
    q=m/p;
 
  /*  fprintf(stderr,"Hello,world! Process %d of %d on %s\n",myid,numprocs,processor_name); */
  
   printf("Hello, Process %d of %d \n",myid,numprocs);
   
   starttime = MPI_Wtime();   //计时开始

   root = 0 ;  //主进程号

   /* 若申请内存不成功，则退出程序 */ 
   if((A = (double *)malloc(doublesize * m * 4))==NULL)
   {
      printf("A can't allocate more memory,terminating. \n");
      exit(1);
   }

   /* 准备数据 */
   n=q; 
   preDate(A,n,myid,numprocs);
   /* 验证数据,打印输出*/
   printMat(A,n,myid,numprocs);
   
   /* 若申请内存不成功，则退出程序 */ 
   if((X = (double *)malloc(doublesize * n))==NULL)
   {
      printf("X can't allocate more memory,terminating. \n");
      exit(1);
   }

   Mat2p = (double *)malloc(doublesize * 2 * p * 4);
   Mat2  = (double *)malloc(doublesize * 2 * 4 );
   X2p   = (double *)malloc(doublesize * 2 * p);
   X2    = (double *)malloc(doublesize * 2 );

/* ########################################################## */

      //分块计算求缩减方程
      n=q;
      Mat=A;
      //===各进程（myid）对子方程进行约化，获取两行的子缩减方程(8个数据)===
      LinearEqRecduce(Mat,n,Mat2);  //对Matq(q,4)对应的三对角方程进行约化
  
      /* 收集从进程发回的数据（8个double值）*/
       root = 0; 
       count = 8;
       MPI_Gather(Mat2,count,MPI_DOUBLE,Mat2p,count,MPI_DOUBLE,root,MPI_COMM_WORLD);
      //求解缩减后的方程组
     if(myid==0)
     { 
      n=2*p;
      forwardLU(Mat2p,n,X2p);  //主进程求解2*p行数的缩减方程，获得2*p个解，x2p(2*p)
      }
     
       /* =========同步 ==========*/
      MPI_Barrier(MPI_COMM_WORLD);

      /* 将缩减方程组的解全发送到从进程（每进程收两个值）*/
       root = 0; 
       count= 2;
       MPI_Scatter(X2p,count,MPI_DOUBLE,X2,count,MPI_DOUBLE,root,MPI_COMM_WORLD);

       /* 利用传回的解进行回代，求出各子块的最终解 */
       n = q;
       Mat = A;
       LinearEqBacktrack(Mat,n,X2,X);   //对Mat(n,4)对应的三对角方程进行回代

/* ########################################################## */

    /* 输出方程的解　*/
    printVal(X,n,myid,numprocs);

   /* 回收内存 */
   recovery(A,X);
   recovery(Mat2p,X2p);
   recovery(Mat2,X2);
   //free(A);                     //回收内存，防止内存溢出
   //free(X);
   

    endtime = MPI_Wtime();   //计时结束
    walltime = endtime - starttime;   //计算总时间
    printf("walltime=%lf\n",walltime);
    
   MPI_Finalize();  /*MPI finalize */
   return 0;

 }
