# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#define PI 3.14


void ndgrid(float* xx, float* yy,float* X, float* Y)
{
	int Res=255,i,j;
	/*for(i=0;i<Res;i++)
	{
		printf("%lf\n",xx[i]);
	}*/
	for(j=0;j<Res;++j)
	{
		for (i=0;i<Res;++i)
		{
			X[(j*Res)+i]=xx[i];
			
		}
	}
	for(i=0;i<Res;++i)
	{
		for (j=0;j<Res;++j)
		{
			Y[(i*Res)+j]=yy[j];
		}
	}
	/*
	for(i=0;i<Res;i++)
	{
		for (j=0;j<Res;j++)
		{
			printf("%lf\t",Y[i][j]);
		}
		printf("\n");
		getchar();
	} 
	
	printf("%lf\n",X[0][0]);
	printf("%lf\n",X[64][64]);
	printf("%lf\n",Y[0][0]);
	printf("%lf\n",Y[64][64]);
	printf("%lf\n",X[24][23]);
	printf("%lf\n",Y[24][23]);
	*/
	//getchar();
	return;
}




