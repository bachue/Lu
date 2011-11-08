###################################################################
#  三对角线性方程组求解的项目文件
#
#  测试主函数：main.c
#  测试函数的头文件和定义文件:test.h,test.c
#  LU求解算法的函数的头文件和定义文件：lu.h lu.c
#  RCD求解算法的函数的头文件和定义文件：rcd.h rcd.c
#  编译与连接：make
#  清除编译与连接的中间文件，只为保留程序的源文件：make　clean
#
##################################################################
objects= test.o lu.o 
#主测试程序
objects+= main.o  
#可运行文件名
execname=trieq_mpi
CC=mpicc

$(execname): $(objects) 
	$(CC) -o $(execname) $(objects) 
trieq:$(objects)
	$(CC) -o trieq $(objects)
main.o:main.c test.h lu.h  trieq.h
	$(CC) -c main.c
test.o:test.c  test.h  trieq.h
	$(CC) -c test.c 
lu.o:lu.c  lu.h  trieq.h
	$(CC) -c lu.c
clean:
	rm *.o trieq_mpi
.PHONY : run
run:
	mpdboot -n 4 -f ~/mpd.hosts
	mpdtrace
	mpirun -np 4 ./$(execname)
	mpdallexit

