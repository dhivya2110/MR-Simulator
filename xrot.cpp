#include <stdio.h>
#include <math.h>

void xrot(float phi,float* Rx)
{
	Rx[0] = 1;
	Rx[1] = 0;
	Rx[2] = 0;
	Rx[3] = 0;
	Rx[4] = cos(phi);
	Rx[5] = -sin(phi);
	Rx[6] = 0;
	Rx[7] = sin(phi);
	Rx[8] = cos(phi);
	/*
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("%lf\t",Rx[i][j]);
		}
		printf("\n");
	}
	getchar();
	*/
	return;
}