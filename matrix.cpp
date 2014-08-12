#include <stdio.h>
#include <stdlib.h>

void matmul(float m1[3][3], float m2[3][3])
{
	int i,j,k;
	float prod[3][3]={0};
	for(i=0;i<3;++i)
	{
    for(j=0;j<3;++j)
		{
	      	for(k=0;k<3;k++)
			{
			prod[i][j]=prod[i][j]+(m1[i][k]*m2[k][j]);
			}
		} 
	}
	for(i=0;i<3;++i)
	{
		for(j=0;j<3;++j)
		{
			m1[i][j]=prod[i][j];
		}
		printf("\n");
	}
	/*
	//Print matrix A
	printf("\nThe value of A is \n");
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("\t%lf\t\t ",prod[i][j]);
		}
		printf("\n");
	}
	*/
	return;
}
void transpose(float* X, int Res)
{
	float* T=(float*)malloc(sizeof(float)*Res*Res);
	int i,j;
	for(i=0;i<Res;++i)
	{
		for (j=0;j<Res;j++)
		{
			T[(i*Res)+j]=X[(j*Res)+i];
		}
	}
	for(i=0;i<Res;++i)
	{
		for (j=0;j<Res;++j)
		{
			X[(i*Res)+j]=T[(i*Res)+j];
			//printf("%lf\t",X[i][j]);
		}
		//printf("\n");
		//getchar();
	}
	free(T);
}