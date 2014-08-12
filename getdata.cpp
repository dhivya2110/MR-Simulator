#include <stdio.h>
#include <stdlib.h>

void getdata(float* c_rhomat, float* c_T1mat, float* c_T2mat)
{
	FILE *rhomat,*t1mat,*t2mat;
	int i,j;
	rhomat =fopen("rhomat_255.dat","r");
	
	/*if(rhomat==NULL)
	{
    printf("Error in Opening file");
    }
	else
	{
		printf("File read successful");
		getchar();
	}*/

	for(i=0;i<255;i++)
	{
		for(j=0;j<255;j++)
		{
			fscanf(rhomat,"%f ",&c_rhomat[(i*255)+j]);
		}
	}
	/*
	for(i=0;i<255;i++)
	{
		for(j=0;j<255;j++)
		{
		printf("\t%f",c_rhomat[(i*255)+j]);
		}
		printf("\n");
	}*/
	fclose(rhomat);
	//getchar();

	t1mat = fopen("t1mat_255.dat","r");
	
	/*if(t1mat==NULL)
	{
    printf("Error in Opening file");
    }
	else
	{
		printf("File read successful");
		getchar();
	}*/

	for(i=0;i<255;i++)
	{
		for(j=0;j<255;j++)
		{
		fscanf(t1mat,"%f ",&c_T1mat[(i*255)+j]);
		}
	}
	/*
	for(i=0;i<255;i++)
	{
		for(j=0;j<255;j++)
		{
		printf("\t%d",c_T1mat[(i*255)+j]);
		}
		printf("\n");
	}*/
	fclose(t1mat);
	//getchar();

	t2mat = fopen("t2mat_255.dat","r");
	
	/*if(t2mat==NULL)
	{
    printf("Error in Opening file");
    }
	else
	{
		printf("File read successful");
		getchar();
	}*/

	for(i=0;i<255;i++)
	{
		for(j=0;j<255;j++)
		{
		fscanf(t2mat,"%f ",&c_T2mat[(i*255)+j]);
		}
	}
	//printf("%lf\n",c_T2mat[(65*255)+65]);
	/*for(i=0;i<255;i++)
	{
		for(j=0;j<255;j++)
		{
		printf("\t%f",c_T2mat[(i*255)+j]);
		}
		printf("\n");
	}*/
	fclose(t2mat);
	//getchar();
	return;
}
