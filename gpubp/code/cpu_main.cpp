#include "common.h"
#include "bp.h"
#include <ctype.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>


#define RUN 1 
#define TEST 0

void reset_network(network* pnet);
void init_network(network* pnet);
void update_node(network* pnet,frame* pframe,int n);
void cpu_bp(network* pnet,int niter, int mode=SUMPRODUCT);

#if RUN
int main(int argc,char* argv[])
{
  network net;
  state input,output;
  reset_network(&net);

  init_messages(&h_src);
  init_messages(&h_dst);

  ifstream fnet("net.txt");
  ifstream fstate("net-input.txt");
  ofstream fout("cpu-output.txt");
  /*get input*/
  fnet>>net;
  cout<<net;

  fstate>>input;
  for(int n=0;n<net.num;n++)
	{
	  if(net.var[n].nchild==0)
		{
		  for(int x=0;x<net.var[n].size;x++)
			net.var[n].lambda[x]=input.value[n][x];
		}
	}
  cout<<input<<endl;

  /*do actual belief prop*/
  int niter=2;
  long start=clock();
  cpu_bp(&net,niter);
  long end=clock();
  printf("CPU time elapsed with %d nodes:%.6f\n",net.num,(float)(end-start)/CLOCKS_PER_SEC);
  free_network_resources(&net);
}
#endif
#if TEST
int main(int argc,char* argv[])
{
  network net;
  //initialize network
  reset_network(&net);
  ifstream fin("out.txt");
  fin>>net;
  cout<<"dumping network"<<endl;
  cout<<net;
  free_network_resources(&net);

  int uidx[3];
  int sizes[]={2,2,2};
  for(int m=0;m<8;m++)
  {
	getidx(sizes,uidx,3,m);
	for(int i=0;i<3;i++)
	  cout<<uidx[i]<<" ";
	cout<<endl;
  }
  fin.close();

  state s;
  ifstream fstate("states.txt");
  fstate>>s;
  cout<<s<<endl;

  network snet;
  reset_network(&snet);
  cout<<"----------------------------------------"<<endl;
  snet.num=3;
  /*sizes*/
  snet.sizes[0]=snet.sizes[1]=snet.sizes[2]=2;
  snet.dag[1][0]=0;
  snet.dag[2][0]=0;
  init_network(&snet);
  init_messages(&h_src);
  init_messages(&h_dst);

  /*CPT*/
  float* M1=(snet.var[1].M);
  float* M2=(snet.var[2].M);
  for(int i=0;i<4;i++)
	M1[i]=M2[i]=0;
  snet.var[1].M[3]=1;
  snet.var[2].M[1]=1;
  /*lambda*/
  float* l1=(snet.var[1].lambda);
  float* l2=(snet.var[2].lambda);
  l1[0]=0;l1[1]=1; l2[0]=1;l2[1]=0;
  for(int iter=0;iter<2;iter++)
	{
	  cout<<"---------------Iteration("<<iter<<")-------------"<<endl;
	  for(int n=0;n<snet.num;n++)
		{
		  update_node(&snet,n);
		  cout<<"node"<<n<<":";
		  for(int i=0;i<2;i++)
			cout<<snet.var[n].bel[i]<<" ";
		  cout<<endl;
		}
	  memcpy(&h_src,&h_dst,sizeof(message));

	}
  free_network_resources(&snet);
}
#endif
