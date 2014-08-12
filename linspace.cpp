# include <stdio.h>
# include <stdlib.h>
# include <math.h>

float* linspace(float min, float max, int n)
{
	float* result=(float*)malloc(sizeof(float)*n);
	int i;
	for (i = 0; i < n; i++)
	{
		float temp = min + i*(max-min)/(floor((float)n) - 1);
		result[i]=temp;
		//printf("%lf\n",result[i]);
	}
	return(result);	
}