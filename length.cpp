#include <stdio.h>
#include <math.h>

int length(float min, float interval, float max)
{
	int Ne=((max-min)/interval);
	if(fmod((max-min),interval)!=0)
	{
		Ne+=1;
	}
	/*
	printf("Min %lf\n",min);
	printf("Interval %lf\n",interval);
	printf("Max %lf\n",max);
	printf("Length %d\n",Ne);
	getchar();
	*/
	return Ne;
}