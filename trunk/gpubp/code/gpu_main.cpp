#include "common.h"
#include "bp.h"
#include <ctype.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "cuda_runtime.h"
#include "cuda.h"
#include "cutil.h"



#define TESTGPU 1
#define TESTCODE 0 

void gpu_bnet(network* pnet_h,int niter,int mode=SUMPRODUCT);
void free_gpu_network_resources(network **ppnet);
void copy_network(network** ppdst,network* psrc,int mode,int alloc=1);


#if TESTCODE
  int main(int argc,char* argv[])
	{
	  ifstream fin("net.txt");
	  network *hnet,*gnet;

	  hnet=new network;
	  fin>>*hnet;
	  cout<<"-----------------INITIAL----------"<<endl;
 	  cout<<*hnet;
	  copy_network(&gnet,hnet,TOGPU,1);
	  free_network_resources(hnet);
	  memset(hnet,0,sizeof(network));
	  cout<<"-----------------FREED----------"<<endl;
 	  cout<<*hnet;
	  memset(hnet,0,sizeof(network));
	  copy_network(&hnet,gnet,FROMGPU,0);
	  cout<<"-----------------COPIED BACK----------"<<endl;
 	  cout<<*hnet;
	 return 0;
	}
#endif

#if TESTGPU
  int main(int argc,char* argv[])
	{
	  ifstream fin("net.txt");
	  ifstream fstate("net-input.txt");

	  network *hnet,*gnet;
	  state   input;

	  hnet=new network;
	  fin>>*hnet;
	  fstate>>input;

	  for(int n=0;n<hnet->num;n++)
	  {
		cout<<"size:"<<hnet->var[n].size<<","<<hnet->sizes[n]<<endl;
		if(hnet->var[n].nchild==0)
		{
		  for(int x=0;x<hnet->var[n].size;x++)
			hnet->var[n].lambda[x]=input.value[n][x];
		}
	 }
	 cout<<*hnet;

	 gpu_bnet(hnet,1);
	 cout<<*hnet;
	 cout<<"------------Beliefs------------"<<endl;
	 for(int n=0;n<hnet->num;n++)
		{
		  cout<<"node"<<n<<":";
		  for(int i=0;i<hnet->var[n].size;i++)
			cout<<hnet->var[n].bel[i]<<" ";
		  cout<<endl;
		}
	 return 0;
	}
#endif
