#include "common.h"
#include "bp.h"
#include <ctype.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "cuda.h"
#include "cutil.h"

#define NTHREADS 4 


__device__ void dgetidx(const int *psizes,int* uidx,int npar,int m)
{
  int stride[MAX_PARENTS];
  stride[0]=1;
  for(int i=1;i<npar;i++)
	stride[i]=stride[i-1]*psizes[i];
  for(int n=(npar-1);n>=0;n--)
	{
	  uidx[n]=m/stride[n];
	  m-=(uidx[n]*stride[n]);
	}
}

void copy_network(network** ppdst,network* psrc,int mode,int alloc=1)
{
  network h_copy;
  if(mode==TOGPU)
	{
	  memcpy(&h_copy,psrc,sizeof(network));

	  for(int i=0;i<psrc->num;i++)
		{
		  if(alloc)
			{
			  float * ptr;
			  cudaMalloc((void**)&(h_copy.var[i].M),psrc->var[i].Msize*sizeof(float));
			}
		    cudaMemcpy(h_copy.var[i].M,psrc->var[i].M,sizeof(float)*psrc->var[i].Msize,cudaMemcpyHostToDevice);
		}
		if(alloc)
		  cudaMalloc((void**)ppdst,sizeof(network));
		cudaMemcpy(*ppdst,&h_copy,sizeof(network),cudaMemcpyHostToDevice);
	}
  else //FROM GPU
	{
	  /*make sure you free resources*/
	  if(alloc)
		*ppdst=(network*)malloc(sizeof(network));
	  cudaMemcpy(*ppdst,psrc,sizeof(network),cudaMemcpyDeviceToHost);
	  
	  network* pdst=*ppdst;
	  for(int i=0;i<pdst->num;i++)
		{
		  float* ptr=pdst->var[i].M;
  		  pdst->var[i].M=(float*)malloc(pdst->var[i].Msize*sizeof(float));
		  cudaMemcpy(pdst->var[i].M,ptr,pdst->var[i].Msize*sizeof(float),cudaMemcpyDeviceToHost);
		}
	}
}

void free_gpu_network_resources(network **ppnet)
{
  network net;
  cudaMemcpy(&net,*ppnet,sizeof(network),cudaMemcpyDeviceToHost);
  for(int i=0;i<net.num;i++)
	{
	  cudaFree(net.var[i].M);
	}
  cudaFree(*ppnet);
  *ppnet=0;
}

typedef frame* frameptr;

__global__ void dummy(frameptr pframe)
{
  for(int i=0;i<MAX_NODES;i++)
  {
    frameptr ptr=&pframe[i];
    ptr->cidx[0]=1;
    ptr->psizes[0]=10;
  }
}
__global__ void update_belief(void* netptr,void* fptr,void* srcptr,void* dstptr)
{
  int tidx=threadIdx.x;
  int n=blockIdx.x;
  /*determine indices from children*/
  int p,c,i,x;
  network* pnet  =(network*)netptr;
  message* g_src =(message*)srcptr;
  message* g_dst =(message*)dstptr;
   frame*  g_frame=(frame*)fptr;

  frame* pframe=&g_frame[n];
  int nchild=pnet->var[n].nchild;
  int npar  =pnet->var[n].npar;
  int size  =pnet->var[n].size;

  /*determine boundaries of each thread*/
  int csize = (nchild+NTHREADS-1)/NTHREADS;
  int cstart= tidx*csize;
  int cstop = min(cstart+csize,nchild);

  int psize = (npar+NTHREADS-1)/NTHREADS;
  int pstart= tidx*psize;
  int pstop = min(pstart+psize,npar);

  int xsize = (size+NTHREADS-1)/NTHREADS;
  int xstart= tidx*xsize;
  int xstop = min(xstart+xsize,size);
  node* np=&(pnet->var[n]);
  
  for(int c=cstart;c<cstop;c++)
	{
	  node* child=&(pnet->var[np->children[c]]);
	  for(p=0;p<child->npar;p++)
	  {
	  	if(child->parents[p]==n)
	    {
		  pframe->cidx[c]=p;
		  break;
		}
	  }
	}
  __syncthreads();
  /*determine indices from parents*/
  for(p=pstart;p<pstop;p++)
	{
	  node* parent=&(pnet->var[np->parents[p]]);
	  for(c=0;c<parent->nchild;c++)
	    {
		  if(parent->children[c]==n)
			{
			  pframe->pidx[p]=c;
			  break;
			}
		}
	}
  __syncthreads();
  /*size of parents*/
  if(tidx==0)
	{
	  for(pframe->msize=1,p=0;p<np->npar;p++)
		{
		  node* parent=&(pnet->var[np->parents[p]]);
		  pframe->psizes[p]=parent->size;
		  pframe->msize*=pframe->psizes[p];
		}
	}
  __syncthreads();

  /*update lambda*/
  for(x=xstart;x<xstop;x++)
	{
	  pframe->lambda[x]=1;
	  if(np->nchild==0)
		{
		  pframe->lambda[x]=np->lambda[x];
		}
	  else
		{
		  for(c=0;c<np->nchild;c++)
			{
			  float cval=g_src->lambda[np->children[c]][pframe->cidx[c]][x];
			  pframe->lambda[x]*=cval;
			}
		}
	}
   __syncthreads();
  /*update pi*/
  for(x=xstart;x<xstop;x++)
	{
	  pframe->pi[x]=0;
	  if(np->npar==0)
		{
		  pframe->pi[x]=np->M[x];
		}
	  else
		{
		  for(int m=0;m<pframe->msize;m++)
			{
			  float prod=1;
			  dgetidx(pframe->psizes,pframe->uidx,np->npar,m);
			  for(p=0;p<np->npar;p++)
				{
				  prod*=g_src->pi[np->parents[p]][pframe->pidx[p]][pframe->uidx[p]];
				}
			  pframe->pi[x]+=np->M[m+x*pframe->msize]*prod;
			}
		}
	}

  __syncthreads();
  float sum=0;

  /*update belief*/
  if(tidx==0)
	{
	  for(x=0;x<np->size;x++)
		{
		  pframe->bel[x]=pframe->lambda[x]*pframe->pi[x];
		  sum+=pframe->bel[x];
		}
	  for(x=0;x<np->size;x++)
		{
		  pnet->var[n].bel[x]=pframe->bel[x]/sum;
	   }
	}
	__syncthreads();
  /*bottom-up propagation*/
  for(p=0;p<np->npar;p++)
	{
	  int   curpar=np->parents[p];
	  node* parent=&(pnet->var[curpar]);
	  sum=0;
	  for(x=xstart;x<xstop;x++)
		{
		  for(int m=0;m<pframe->msize;m++)
			{
			  float prod=1;
			  dgetidx(pframe->psizes,pframe->uidx,np->npar,m);
			  for(int pp=0;pp<np->npar;pp++)
				{
				  int thispar=np->parents[pp];
				  if(pp==p)continue;
				  prod*=g_src->pi[thispar][pframe->pidx[pp]][pframe->uidx[pp]];
				}
			  pframe->term[n][pframe->uidx[p]][x]+=np->M[x*pframe->msize+m]*prod;
			}
		}
	  __syncthreads();/*wait for everybody to compute sum*/
	  float usum=0;
	  int u=0;
	  if(tidx==0)
		{
		  for(u=0;u<parent->size;u++)
			{
			  sum=0;
			  for(x=0;x<np->size;x++)
				{
				  sum+=pframe->term[n][u][x]*pframe->lambda[x];
				}
			  usum+=sum;
			  g_dst->lambda[n][p][u]=sum;
			}
		  for(int u=0;u<parent->size;u++)
			{
			  g_dst->lambda[n][p][u]/=usum;
			}
		}
	}/*end parent*/
    __syncthreads();
   //top down propagation
	for(c=cstart;c<cstop;c++)
	  {
		int curchild=np->children[c];
		sum=0;
		for(x=0;x<np->size;x++)
		  {
			pframe->py[x]=pframe->bel[x]/(g_src->lambda[curchild][pframe->cidx[c]][x]+1e-5);
			sum+=pframe->py[x];
		  }
		for(x=0;x<np->size;x++)
		  {
			g_dst->pi[n][c][x]=pframe->py[x]/sum;
		  }
	  }
	__syncthreads();
}


void gpu_bnet(network* pnet_h,int niter,int mode=SUMPRODUCT)
{
  network* pnet_g; //gpu_copy
  message* g_src,*g_dst; /*gpu messages*/
  frame *g_frame; /*gpu frames*/
  int f,n,x,c,p,m;

  cout<<"Initializing messages"<<endl;
  /*allocate message structures*/
  CUDA_SAFE_CALL(cudaMalloc((void**)&g_src,sizeof(message)));
  CUDA_SAFE_CALL(cudaMalloc((void**)&g_dst,sizeof(message)));

  init_messages(&h_src); /*initialize on host*/
  CUDA_SAFE_CALL(cudaMemcpy(g_src,&h_src,sizeof(message),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(g_dst,&h_src,sizeof(message),cudaMemcpyHostToDevice));
  cout<<"done"<<endl;

  /*allocate_frame*/
  CUDA_SAFE_CALL(cudaMalloc((void**)&g_frame,MAX_NODES*sizeof(frame)));
  CUDA_SAFE_CALL(cudaMemset(g_frame,0,MAX_NODES*sizeof(frame)));
  /*copy to gpu*/
  copy_network(&pnet_g,pnet_h,TOGPU,1);
  long start_time=clock();
  for(int iter=0;iter<2;iter++)
	{
	  /*update belief*/
	  dim3 gridSize(3,1,1);
	  dim3 blkSize(10,1,1);
      //dummy<<<gridSize,blkSize>>>(g_frame);
      update_belief<<<gridSize,blkSize>>>(pnet_g,g_frame,g_src,g_dst);
	  CUT_CHECK_ERROR("something bad happened");
	  /*update messages*/
	  CUDA_SAFE_CALL(cudaThreadSynchronize());
	  CUDA_SAFE_CALL(cudaMemcpy(g_dst,g_src,sizeof(message),cudaMemcpyDeviceToDevice));
  #if 1
      cout<<"--------Frames("<<iter<<")---------"<<endl;
      cudaMemcpy(h_frame,g_frame,sizeof(frame)*MAX_NODES,cudaMemcpyDeviceToHost);

    for(n=0;n<pnet_h->num;n++)
    { 
	    cout<<"-----------------------------------node("<<n<<")"<<endl;
	    int  sz=pnet_h->var[n].size;

	    cout<<"size:"<<sz<<endl;

        cout<<"cidx:";
	    for(c=0;c<pnet_h->var[n].nchild;c++)
		    cout<<h_frame[n].cidx[c]<<" ";
    	cout<<endl;

	    cout<<"pidx:";
	    for(p=0;p<pnet_h->var[n].npar;p++)
	    	cout<<h_frame[n].pidx[p]<<" ";
	    cout<<endl;

    	cout<<"psize:";
	    for(p=0;p<pnet_h->var[n].npar;p++)
		    cout<<h_frame[n].psizes[p]<<" ";
    	cout<<endl;

	    cout<<"msize:"<<h_frame[n].msize<<endl;

	    cout<<"lambda:";
    	for(x=0;x<sz;x++)
	    {
		    cout<<h_frame[n].lambda[x]<<" ";
	    }
	    cout<<endl;	
  	    cout<<"pi:";
	    for(x=0;x<sz;x++)
	    {
		    cout<<h_frame[n].pi[x]<<" ";
	    }
	    cout<<endl;	
    }	   
    cout<<"---------Beliefs-------"<<endl;
    for(n=0;n<pnet_h->num;n++)
    {
	  cout<<"node"<<n<<":";
	  for(int i=0;i<pnet_h->var[n].size;i++)
		{
		  cout<<pnet_h->var[n].bel[i]<<" ";
		}
	  cout<<endl;
    }
        #endif
  }
  long end_time=clock();
  printf("GPU time elapsed with %d nodes:%.6f\n",pnet_h->num,(float)(end_time-start_time)/CLOCKS_PER_SEC);

  /*copy back to cpu*/
  free_network_resources(pnet_h);
  memset(pnet_h,0,sizeof(network));
  copy_network(&pnet_h,pnet_g,FROMGPU,0);
  cout<<*pnet_h;

 /*clean up*/
  free_gpu_network_resources(&pnet_g);
  cudaFree(g_frame);
  cudaFree(g_src);
  cudaFree(g_dst);
}

