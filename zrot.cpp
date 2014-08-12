#include <math.h>

void zrot(float phi,float Rz[3][3])
{
	Rz[0][0] = cos(phi);
	Rz[0][1] = -sin(phi);
	Rz[0][2] = 0;
	Rz[1][0] = sin(phi);
	Rz[1][1] = cos(phi) ;
	Rz[1][2] = 0;
	Rz[2][0] = 0;
	Rz[2][1] = 0;
	Rz[2][2] = 1;
	/*for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			printf("%lf\t",Ry[i][j]);
		}
		printf("\n");
	}
	getchar();*/
	return;
}