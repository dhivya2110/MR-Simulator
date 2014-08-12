#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "length.h"
#include "linspace.h"
#include "ndgrid.h"
#include "matrix.h"
#include "getdata.h"
#include "xrot.h"
#include "yrot.h"

#define Res 255
#define PI 3.14
#define Nf 1

// for cuda error checking
#define cudaCheckErrors(msg) \
	do { \
		cudaError_t __err = cudaGetLastError(); \
		if (__err != cudaSuccess) { \
			fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
					msg, cudaGetErrorString(__err), \
					__FILE__, __LINE__); \
			fprintf(stderr, "*** FAILED - ABORTING\n"); \
		} \
	} while (0)

typedef struct complex{
	float real;
	float imag;
}complex;

void sliceout_spinechoV2(complex* K0,int lam, int Ne, int NRep,float* h_X, float* h_Y, float* h_rhomat, float* h_t1mat, float* h_t2mat, float* h_Rflipx,float* h_Rflipy,float Gx,float Gy,float Gama, float dT, float TRep);

__constant__ float d_Gx;
__constant__ float d_Gy;
__constant__ int d_Ne;
__constant__ int d_NRep;
__constant__ int d_lam;
__constant__ float d_Gama;
__constant__ float d_dT;
__constant__ float d_TRep;

__device__ void free_precess(float dTf, float T1, float T2, float dF, float* A, float* B)
{
	//Declare local variables
	int i,j,k;
	float E1,E2,phi;
	//Calculate phi, E1 and E2
	phi=2*PI*dF*dTf/1000;
	E1=exp(-dTf/T1);
	E2=exp(-dTf/T2);
	//Allot values to matrix A and B
	A[0] =E2;
	A[1] =0;	
	A[2] =0;
	A[3] =0;
	A[4] =E2;
	A[5] =0;
	A[6] =0;
	A[7] =0;
	A[8] =E1;

	B[0] =0;
	B[1] =0;
	B[2] =1-E1;

	float Rz[3][3]={0};

	Rz[0][0] = cos(phi);
	Rz[0][1] = -sin(phi);
	Rz[0][2] = 0;
	Rz[1][0] = sin(phi);
	Rz[1][1] = cos(phi) ;
	Rz[1][2] = 0;
	Rz[2][0] = 0;
	Rz[2][1] = 0;
	Rz[2][2] = 1;
	float prod[3][3]={0};
	for(i=0;i<3;++i)
	{
		for(j=0;j<3;++j)
		{
			for(k=0;k<3;k++)
			{
				prod[i][j]=prod[i][j]+(A[(i*3)+k]*Rz[k][j]);
			}
		} 
	}
	for(i=0;i<3;++i)
	{
		for(j=0;j<3;++j)
		{
			A[(i*3)+j]=prod[i][j];
		}
	}
	return;
}

__device__ void spinecho(float x, float y, float T1, float T2, float* Morg, float* d_Mout, float* d_Rflipx, float* d_Rflipy, complex* d_K0, int id)
{

	//printf("Gy = %f in device\n", d_Gy);
	int i,j,k,s=d_Ne-(d_lam+1),s1=d_Ne-d_lam-1,s2=d_Ne-1,s3=(3*d_Ne)-1;
	float dF=0,temp=0,M0[3]={0,0,0},ans[3]={0,0,0};
	//cuPrintf("Gama %f",d_Gama);
	//cuPrintf("Gx %f", d_Gx);
	//cuPrintf("Gy %f",d_Gy);
	//cuPrintf("t2 %f \n",T2);
	//printf("morg[3] = %f", Morg[2]);
	float dFx=(dF+(d_Gama*x*d_Gx)/(2*PI));
	float dFy=(dF+(d_Gama*y*d_Gy)/(2*PI));
	float dTf=(d_TRep-(((3*d_Ne)+1)*d_dT));
	float A0[9],B0[3],Ax[9],Bx[3],Ay[9],By[3],Af[9],Bf[3];
	float M[3*766];

	free_precess(d_dT,T1,T2,dF,A0,B0);
	free_precess(d_dT,T1,T2,dFx,Ax,Bx);
	free_precess(d_dT,T1,T2,dFy,Ay,By);
	free_precess(dTf,T1,T2,dF,Af,Bf);
	/*for(i=0;i<3;i++)
	  {
	  d_Mout[i]=d_Rflipy[i];
	  }*/
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			temp=temp+d_Rflipy[(i*3)+j]*Morg[j];
		}
		ans[i]=temp;
		temp=0;
		M[(i*d_NRep)+0]=ans[i];
		ans[i]=0;
		//d_Mout[i]=M[(i*d_NRep)+0];
	}

	for(k=1;k<s;k++)
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				temp=temp+Ax[(i*3)+j]*M[(j*d_NRep)+(k-1)];
			}
			ans[i]=temp;
			temp=0;
			M[(i*d_NRep)+k]=ans[i]+Bx[i];
			ans[i]=0;
			//d_Mout[i]=M[(i*d_NRep)+k];
		}
	}

	for(k=s1;k<s2;k++)
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				temp=temp+Ay[(i*3)+j]*M[(j*d_NRep)+(k-1)];
			}
			ans[i]=temp;
			temp=0;
			M[(i*d_NRep)+k]=ans[i]+By[i];
			//d_Mout[i]=M[(i*d_NRep)+k];
		}
	}

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			temp=temp+d_Rflipx[(i*3)+j]*M[(j*d_NRep)+(d_Ne-2)];
		}
		M[(i*d_NRep)+(d_Ne-1)]=temp;
		temp=0;
		//d_Mout[i]=M[(i*d_NRep)+(d_Ne-1)];
	}

	for(k=d_Ne;k<(3*d_Ne);k++)
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				temp=temp+Ax[(i*3)+j]*M[(j*d_NRep)+(k-1)];
			}
			ans[i]=temp;
			temp=0;
			M[(i*d_NRep)+k]=ans[i]+Bx[i];
			ans[i]=0;
			//d_Mout[i]=M[(i*d_NRep)+k];
		}
	}

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			temp=temp+Af[(i*3)+j]*M[(j*d_NRep)+s3];
		}
		ans[i]=temp;
		temp=0;
		M0[i]=ans[i]+Bf[i];
		ans[i]=0;
	}

	for(i=0;i<3;i++)
	{
		M0[i]=(1.0/Nf)*M0[i];
		//d_Mout[i]=M0[i];
	}

	for(i=0;i<3;i++)
	{
		for(j=0;j<d_NRep;j++)
		{
			M[(i*d_NRep)+j]=(1.0/Nf)*M[(i*d_NRep)+j];
			//cuPrintf("M Real in loop %f \n",M[(i*d_NRep)+j]);
		}
	}

	__syncthreads();

	for(i=d_Ne;i<(3*d_Ne);i++)
	{
		atomicAdd(&d_K0[i-d_Ne].real,M[i]);
		atomicAdd(&d_K0[i-d_Ne].imag,M[i+d_NRep]);
	}
	return;
}

__global__ void spinKernel(float* d_rhomat, float* d_t1mat, float* d_t2mat, float* d_X,float* d_Y,float* d_Rflipx,float* d_Rflipy,float* d_Mout,complex* d_K0)
{
	//printf("ne %d\n",d_Ne);
	//printf("nrep %d\n",d_NRep);
	int i;
	//int i=blockIdx.x;
	int j=threadIdx.x;
	int tid=(blockIdx.x*blockDim.x)+threadIdx.x;
	float x=d_X[tid];
	float y=d_Y[tid];
	float Morg[3]={0,0,d_rhomat[tid]};
	float T1=d_t1mat[tid];
	float T2=d_t2mat[tid];
	//printf("Gy = %f in global func\n", d_Gy);
	//float x=1;
	//float y=1;
	//float Morg[3]={1,0,0};
	//float T1= 100;
	//float T2=80;
	spinecho(x,y,T1,T2,Morg,d_Mout,d_Rflipx,d_Rflipy,d_K0,threadIdx.x);
}

int main()
{
	int Fov[2]={-10,10},y0=Fov[1],i,j,Nx=0,N_ky=0,lam=1;
	float flip=90,TE=30.0,dT=TE/Res,Gama=(47*pow(10.0,6)),TRep=(3*TE),Gx=2*PI*pow(10.0,3)/(Gama*dT*(Fov[1]-Fov[0])),delgy=0,Gy=0,max=0;
	int Ne=length(0,dT,(TE/2)),NRep=length(0,dT,TRep);
	float *xx=linspace(Fov[0],Fov[1],Res);
	float *yy=linspace(y0,-y0,Res);

	float* h_X=(float*)malloc(sizeof(float)*Res*Res);
	float* h_Y=(float*)malloc(sizeof(float)*Res*Res);
	ndgrid(xx,yy,h_X,h_Y);
	transpose(h_X,Res);
	transpose(h_Y,Res);

	float* h_rhomat=(float*)malloc(sizeof(float)*Res*Res);
	float* h_t1mat=(float*)malloc(sizeof(float)*Res*Res);
	float* h_t2mat=(float*)malloc(sizeof(float)*Res*Res);
	getdata(h_rhomat,h_t1mat,h_t2mat);

	float* h_Rflipy=(float*)malloc(3*3*sizeof(float));
	float* h_Rflipx=(float*)malloc(3*3*sizeof(float));
	yrot(flip,h_Rflipy);
	xrot(PI,h_Rflipx);

	//printf("Ne %d\n",Ne);
	//printf("NRep %d\n",NRep);

	complex* K0=(complex*)malloc(256*sizeof(complex));
	
	sliceout_spinechoV2(K0,lam,Ne,NRep,h_X,h_Y,h_rhomat,h_t1mat,h_t2mat,h_Rflipx,h_Rflipy,Gx,Gy,Gama,dT,TRep);
	cudaCheckErrors("sliceout_spinecho failed!");
	FILE* output=fopen("sliceout_0_255_cuda.dat","wb");
	for(i=0;i<255;++i)
	{
		fprintf(output,"%f %f\n ",K0[i].real,K0[i].imag);	
	}
	fclose(output);
	printf("File written success\n");

	float* val=(float*)malloc(256*sizeof(float));
	for(i=0;i<255;i++)
	{
		val[i]=sqrt((pow(K0[i].real,2))+(pow(K0[i].imag,2)));
		//printf("%lf\n",val[i]);
	}
	//max of abs
	max=val[0];
	for(i=0;i<255;i++)
	{
		if(max<val[i])
		{
			max=val[i];
			Nx=i+1;
		}
	}
	free(val);
	//free(K0);
	N_ky=(2*Nx)-1;
	delgy=Gx/lam;
	Gy=(Nx-1)*delgy;
	printf("Nx = %d\n",Nx);
	printf("N_ky = %d\n", N_ky);
	//getchar();
	float* gy=linspace(Gy,-Gy,(2*Nx)-1);
	complex* K=(complex*)malloc(N_ky*N_ky*sizeof(complex));
	for(j=0;j<N_ky;j++)
	{
		sliceout_spinechoV2(K0,lam,Ne,NRep,h_X,h_Y,h_rhomat,h_t1mat,h_t2mat,h_Rflipx,h_Rflipy,Gx,gy[j],Gama,dT,TRep);
		cudaCheckErrors("sliceout_spinecho failed!");
		//K(j,:)=k(1: N_ky);
		for(i=0;i<N_ky;i++)
		{
			K[(j*N_ky)+i]=K0[i];
		}
		
	/*	char fname[100];
		snprintf(fname, 88888888, "file_%d.txt", j);
		FILE* output=fopen(fname,"wb");
	        for(i=0;i<255;++i)
        	{
                	fprintf(output,"%f %f\n ",K0[i].real,K0[i].imag);
        	}
		fclose(output);*/

		printf("Slice %d completed\n",j);
	}

	FILE* out=fopen("space_255_cuda.dat","wb");
	for(i=0;i<N_ky;i++)
	{
		for(j=0;j<N_ky;j++)
		{
			printf("writing K[%d]", (i*N_ky)+j);
			fprintf(out,"%f %f \n",K[(i*N_ky)+j].real,K[(i*N_ky)+j].imag);	

		}
	}

	fclose(out);
	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaDeviceReset();
	cudaCheckErrors("Device Reset failed!");

	//free(K);
	//free(gy);
	free(K0);
	free(h_Rflipx);
	free(h_Rflipy);
	free(h_t2mat);
	free(h_t1mat);
	free(h_rhomat);
	free(h_Y);
	free(h_X);
	free(yy);
	free(xx);
	return 0;
}

//Helper function to implement CUDA kernel
void sliceout_spinechoV2(complex* K0, int lam, int Ne, int NRep,float* h_X, float* h_Y,float* h_rhomat, float* h_t1mat, float* h_t2mat, float* h_Rflipx, float* h_Rflipy,float Gx,float Gy,float Gama, float dT, float TRep)
{
	float Mout[3];
	float* d_X=0;
	float* d_Y=0;
	float* d_rhomat=0;
	float* d_t1mat=0;
	float* d_t2mat=0;
	float* d_Rflipx=0;
	float* d_Rflipy=0;
	float* d_Mout=0;
	complex* d_K0=0;
	size_t size=Res*Res;
	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaSetDevice(0);
	cudaCheckErrors("Set device failed!");

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaMalloc((void**)&d_X, size*sizeof(float));
	cudaCheckErrors("malloc d_X failed!");

	cudaMalloc((void**)&d_Y, size * sizeof(float));
	cudaCheckErrors("malloc d_Y failed!");

	cudaMalloc((void**)&d_rhomat,size * sizeof(float));
	cudaCheckErrors("malloc rhomat failed!");

	cudaMalloc((void**)&d_t1mat,size * sizeof(float));
	cudaCheckErrors("malloc t1mat failed!");

	cudaMalloc((void**)&d_t2mat,size * sizeof(float));
	cudaCheckErrors("malloc t2mat failed!");

	cudaMalloc((void**)&d_Rflipx, 9* sizeof(float));
	cudaCheckErrors("malloc d_Rflipx failed!");

	cudaMalloc((void**)&d_Rflipy, 9* sizeof(float));
	cudaCheckErrors("malloc d_Rflipy failed!");

	cudaMalloc((void**)&d_Mout, 3* sizeof(float));
	cudaCheckErrors("malloc d_Mout failed!");

	cudaMalloc((void**)&d_K0, 256* sizeof(complex));
	cudaCheckErrors("malloc d_K0 failed!");

	// Copy input vectors from host memory to GPU buffers.
	cudaMemcpy(d_X, h_X, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy X failed!");

	cudaMemcpy(d_Y, h_Y, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy Y failed!");

	cudaMemcpy(d_rhomat, h_rhomat, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy rhomat failed!");

	cudaMemcpy(d_t1mat, h_t1mat, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy rhomat failed!");

	cudaMemcpy(d_t2mat, h_t2mat, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy rhomat failed!");

	cudaMemcpy(d_Rflipx, h_Rflipx, 9 * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy Rflipx failed!");

	cudaMemcpy(d_Rflipy, h_Rflipy, 9 * sizeof(float), cudaMemcpyHostToDevice);
	cudaCheckErrors("memcpy Rflipy failed!");

	cudaMemset(d_Mout,0,3*sizeof(float));
	cudaCheckErrors("memcpy Mout failed!");

	cudaMemset(d_K0,0,256*sizeof(complex));
	cudaCheckErrors("memcpy K0 failed!");

	cudaMemcpyToSymbol(d_Gx, &Gx, sizeof(float));
	cudaMemcpyToSymbol(d_Gy, &Gy, sizeof(float));
	cudaMemcpyToSymbol(d_Ne, &Ne, sizeof(int));
	cudaMemcpyToSymbol(d_NRep, &NRep, sizeof(int));
	cudaMemcpyToSymbol(d_lam, &lam, sizeof(int));
	cudaMemcpyToSymbol(d_Gama, &Gama, sizeof(float));
	cudaMemcpyToSymbol(d_dT, &dT, sizeof(float));
	cudaMemcpyToSymbol(d_TRep, &TRep, sizeof(float));

	// Launch a kernel on the GPU with one thread for each element.
	spinKernel<<<255,255>>>(d_rhomat,d_t1mat,d_t2mat,d_X,d_Y,d_Rflipx,d_Rflipy,d_Mout,d_K0);

	cudaCheckErrors("Kernel launch failed!");

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaDeviceSynchronize();
	cudaCheckErrors("Device Synchronize failed!");

	// Copy output vector from GPU buffer to host memory.
	cudaMemcpy(K0, d_K0, 256 * sizeof(complex), cudaMemcpyDeviceToHost);
	cudaCheckErrors("Memcpy output failed!");

	cudaMemcpy(Mout, d_Mout, 3 * sizeof(float), cudaMemcpyDeviceToHost);
	cudaCheckErrors("Memcpy Mout failed!");

	//printf(" %f %f %f \n ", Mout[0], Mout[1], Mout[2]);
	cudaFree(d_K0);
	cudaFree(d_Mout);
	cudaFree(d_Rflipy);
	cudaFree(d_Rflipx);
	cudaFree(d_t2mat);
	cudaFree(d_t1mat);
	cudaFree(d_rhomat);
	cudaFree(d_Y);
	cudaFree(d_X);
	return;
}
