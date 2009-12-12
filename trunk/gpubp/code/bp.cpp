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


using namespace std;
message h_src,h_dst; /*host messages*/
frame h_frame[MAX_NODES]; /*host frames*/


/*
  format
  node# size val0 val1...
*/
ostream& operator<<(ostream& out,const state& st)
{
  int num=st.num;int i;
  out<<num<<endl;
  for(int i=0;i<num;i++)
	{
	  out<<(i+1)<<" "<<st.sizes[i]<<" ";
	  for(int j=0;j<st.sizes[i];j++)
		{
		  out<<st.value[i][j]<<" ";
		}
	  out<<endl;
	}
  out<<endl;
}

istream& operator>>(istream& in,state& st)
{
  int num;int size;int nstates;
  memset(&st,0,sizeof(state));
  in>>nstates;
  st.num=nstates;
  for(int i=0;i<nstates;i++)
	{
	  in>>num>>size;
	  st.sizes[i]=size;
	  for(int j=0;j<size;j++)
		in>>st.value[i][j];
	}
}


void skipcomments(istream &in)
{
  char buff[1024]={0};
  while(1)
  {
	in>>ws; /*eats white space*/
	if(in.peek()=='#')
	  {
		in.getline(buff,sizeof(buff));
		//cout<<buff<<endl;
	  }
	else
	  break;
  }
}

ostream & operator<<(ostream & out,const node& nd)
{
  out<<"#size"<<endl;
  out<<nd.size<<endl;

  out<<"#parents"<<endl;
  out<<nd.npar<<" ";
  for(int p=0;p<nd.npar;p++)
	out<<(nd.parents[p]+1)<<" ";
  out<<endl;

  out<<"#children"<<endl;
  out<<nd.nchild<<" ";
  for(int c=0;c<nd.nchild;c++)
	out<<(nd.children[c]+1)<<" ";
  out<<endl;
  out<<"#M"<<endl;
  out<<nd.Msize<<" ";
  if(nd.M)
	{
	  for(int m=0;m<nd.Msize;m++)
		out<<setprecision(6)<<nd.M[m]<<" ";
	  out<<endl;
	}
  out<<endl;
}

ostream & operator<<(ostream & out,const network& net)
{
  out<<"#num nodes"<<endl;
  out<<net.num;
  for(int i=0;i<net.num;i++)
  {
	out<<"#----------------(node"<<(i+1)<<")"<<endl;
	out<<net.var[i];
  }
}

istream & operator>>(istream& in,node& node)
{
  skipcomments(in);
  in>>node.size;

  /*initialize node values*/
  for(int x=0;x<node.size;x++)
	node.lambda[x]=node.pi[x]=node.bel[x]=1.00;

  skipcomments(in);
  in>>node.npar;


  skipcomments(in);
  for(int p=0;p<node.npar;p++)
	{
	in>>node.parents[p];node.parents[p]--;
	}
  
  skipcomments(in);
  in>>node.nchild;
  for(int c=0;c<node.nchild;c++)
	{
	in>>node.children[c];node.children[c]--;
	}
  int msize=0;
  float tmp=0;
  skipcomments(in);
  in>>node.Msize;
  if(!node.M)
  {
	node.M=(float*)malloc(node.Msize*sizeof(float));
  }
  for(int m=0;m<node.Msize;m++)
	  in>>node.M[m];
  return in;
}

istream & operator>>(istream& in,network& net)
{
  skipcomments(in);
  in>>net.num;
  cout<<net.num<<endl;
  for(int n=0;n<net.num;n++)
  {
	in>>net.var[n];
	net.sizes[n]=net.var[n].size;
  }
  return in;
}


void reset_network(network* pnet)
{
  memset(pnet->var,0,sizeof(node)*MAX_NODES);
  for(int s=0;s<MAX_NODES;s++)
	pnet->sizes[s]=-1;
  for(int i=0;i<MAX_NODES;i++)
	{
	  pnet->var[i].M=0;
	  for(int j=0;j<MAX_PARENTS;j++)
		{
		  pnet->dag[i][j]=-1;
		}
	}
}

void free_network_resources(network* pnet)
{
  for(int i=0;i<pnet->num;i++)
	{
	  if(pnet->var[i].M)
		{
		  free(pnet->var[i].M);
		  pnet->var[i].Msize=0;
		  pnet->var[i].M=0;
		}
	}
}

/*!
  @function init_network
  @desc initialized the node variables based on dag and sizes 
        assumes network is zeroed out earlier. 
 */
void init_network(network* pnet)
{
  /*initialize */
  int x;
  for(int n=0;n<pnet->num;n++)
	{
	  int msize=1;
	  pnet->var[n].size=pnet->sizes[n];
	  for(int i=0;pnet->dag[n][i]>=0;i++)
		{
		  int parent=pnet->dag[n][i];
		  msize*=pnet->sizes[parent];
		  pnet->var[n].parents[pnet->var[n].npar++]=parent;
		  pnet->var[parent].children[pnet->var[parent].nchild++]=n;
		}
	  pnet->var[n].observed=0;
	  pnet->var[n].Msize=msize*pnet->var[n].size;
	  /*set up lambda,pi,belief*/
	  for(x=0;x<pnet->var[n].size;x++)
      {
		  pnet->var[n].lambda[x]=
		  pnet->var[n].pi[x]=
          pnet->var[n].bel[x]=1.00;
	  }

	  /*set up CPT*/
	  pnet->var[n].M=(float*)malloc(sizeof(float)*pnet->var[n].Msize);
	  memset(pnet->var[n].M,0,sizeof(float)*pnet->var[n].Msize);
	  for(x=0;x<pnet->var[n].size;x++)
		{
		  float val=1.00/pnet->var[n].size;
		  pnet->var[n].lambda[x]=pnet->var[n].pi[x]=val;
		  for(int i=0;i<msize;i++)
			pnet->var[n].M[x*msize+i]=val;
		}
	}
}

void init_messages(message* pm)
{
  	  unsigned long lsize=MAX_NODES*MAX_PARENTS*MAX_STATES;
	  fill((float*)&pm->lambda,(float*)(pm->lambda)+lsize,1.00);
	  unsigned long psize=MAX_NODES*MAX_CHILDREN*MAX_STATES;
	  fill((float*)pm->pi,(float*)(pm->pi)+psize,1.00);
}

void getidx(const int *psizes,int* uidx,int npar,int m)
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


void update_node(network* pnet,frame* pframe,int n)
{
  int msize;

  node* np=&(pnet->var[n]);
  /*determine indices from children*/
  int p,c,i,x;

  for(c=0;c<np->nchild;c++)
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
  /*determine indices from parents*/
  for(p=0;p<np->npar;p++)
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
  /*size of parents*/
  for(pframe->msize=1,p=0;p<np->npar;p++)
	{
	  pframe->psizes[p]=pnet->var[np->parents[p]].size;
	  pframe->msize*=pframe->psizes[p];
	}

  /*update lambda*/
  cout<<"lambda("<<n<<"):";
  for(x=0;x<np->size;x++)
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
			  float cval=h_src.lambda[np->children[c]][pframe->cidx[c]][x];
			  cout<<"("<<cval<<"),";
			  pframe->lambda[x]*=cval;
			}
		}
	  cout<<pframe->lambda[x]<<" ";
	}
  cout<<endl;

  /*update pi*/
  cout<<"pi("<<n<<"):";
  for(x=0;x<np->size;x++)
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
			  getidx(pframe->psizes,pframe->uidx,np->npar,m);
			  for(p=0;p<np->npar;p++)
				{
				  prod*=h_src.pi[np->parents[p]][pframe->pidx[p]][pframe->uidx[p]];
				}
			  pframe->pi[x]+=np->M[m+x*pframe->msize]*prod;
			}
		}
	  cout<<pframe->pi[x]<<" ";
	}
  cout<<endl;
  float sum=0;

  cout<<"bel("<<n<<"):";
  /*update belief*/
  for(x=0;x<np->size;x++)
	{
	  pframe->bel[x]=pframe->lambda[x]*pframe->pi[x];
	  sum+=pframe->bel[x];
	}
  for(x=0;x<np->size;x++)
	{
	  pnet->var[n].bel[x]=pframe->bel[x]/sum;
	  cout<<pframe->bel[x]<<" ";
	}
  cout<<endl;

  /*bottom-up propagation*/
  cout<<"bottom up("<<n<<"):"<<endl;
  for(p=0;p<np->npar;p++)
	{
	  cout<<"to parent("<<np->parents[p]<<"):";
	  int   curpar=np->parents[p];
	  node* parent=&(pnet->var[curpar]);
	  sum=0;
	  for(x=0;x<np->size;x++)
		{
		  for(int m=0;m<pframe->msize;m++)
			{
			  float prod=1;
			  getidx(pframe->psizes,pframe->uidx,np->npar,m);
			  for(int pp=0;pp<np->npar;pp++)
				{
				  int thispar=np->parents[pp];
				  if(pp==p)continue;
				  prod*=h_src.pi[thispar][pframe->pidx[pp]][pframe->uidx[pp]];
				}
			  pframe->term[n][pframe->uidx[p]][x]+=np->M[x*pframe->msize+m]*prod;
			}
		}
	  float usum=0;
	  int u=0;
	  for(u=0;u<parent->size;u++)
		{
		  sum=0;
		  for(x=0;x<np->size;x++)
			{
			  sum+=pframe->term[n][u][x]*pframe->lambda[x];
			}
		  usum+=sum;
		  h_dst.lambda[n][p][u]=sum;
		}
	  for(int u=0;u<parent->size;u++)
		{
  		  h_dst.lambda[n][p][u]/=usum;
		  cout<< h_dst.lambda[n][p][u]<<" ";
		}
	  cout<<endl;
	}

  cout<<"top-down("<<n<<"):"<<endl;
  /*top down propagation*/
  for(c=0;c<np->nchild;c++)
  {
	cout<<"to child("<<np->children[c]<<"):";
	int curchild=np->children[c];
	sum=0;
	for(x=0;x<np->size;x++)
	  {
		pframe->py[x]=pframe->bel[x]/(h_src.lambda[curchild][pframe->cidx[c]][x]+1e-5);
		sum+=pframe->py[x];
	  }
	for(x=0;x<np->size;x++)
	  {
		h_dst.pi[n][c][x]=pframe->py[x]/sum;
		cout<<h_dst.pi[n][c][x]<<" ";
	  }
	cout<<endl;
  }
  cout<<endl<<"----------------------------------"<<endl;
}

/*!
  @function cpu_bp
  @desc     runs belief propagation
  @param    pnet [IN|OUT] 
  @param    model [IN]    chooses between max prod
*/
void cpu_bp(network* pnet,int niter, int mode=SUMPRODUCT)
{
  for(int iter=0;iter<niter;iter++)
	{
	  cout<<"---------------Iteration("<<iter<<")-------------"<<endl;
	  for(int n=0;n<pnet->num;n++)
		{
		  update_node(pnet,&h_frame[n],n);
		  cout<<"node"<<n<<":";
		  for(int i=0;i<pnet->var[n].size;i++)
			cout<<pnet->var[n].bel[i]<<" ";
		  cout<<endl;
		}
	  memcpy(&h_src,&h_dst,sizeof(message));
	}
}

