#include "lu.h"

int forwordLu(float *A, float *X, int N)
{
	int i, j;
	float t;

	c(1) = c(1) / b(1);
	d(1) = d(1) / b(1);

	for(i = 2; i <= (N - 1); i++)
	{
		t = b(i) - a(i) * c(i - 1);
		c(i) = c(i) / t;
		d(i) = (d(i) - (a(i) * d(i - 1))) / t;
	}

	i = N;
	t = b(i) - a(i) * c(i - 1);
	d(i) = (d(i) - (a(i) * d(i - 1))) / t;

	x(N) = d(N);
	for(i = N - 1; i >= 1; i--)
	{
		x(i) = d(i) - c(i) * x(i + 1);
	}
	
	return 0;
}
