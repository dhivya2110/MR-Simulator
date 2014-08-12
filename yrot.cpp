#include <stdio.h>
#include <math.h>

void yrot(float phi,float* Ry)
{
	Ry[0] = cos(phi);
	Ry[1] = 0;
	Ry[2] = sin(phi);
	Ry[3] = 0;
	Ry[4] = 1;
	Ry[5] = 0;
	Ry[6] = -sin(phi);
	Ry[7] = 0;
	Ry[8] = cos(phi);
	/*
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("%lf\t",Ry[i][j]);
		}
		printf("\n");
	}
	getchar();
	*/
	return;
}