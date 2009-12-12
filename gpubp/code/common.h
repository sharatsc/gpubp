#ifndef __COMMON_H__
#define __COMMON_H__

#include <cuda.h>
#include <host_defines.h>

#define MAX_PARENTS  8
#define MAX_CHILDREN 8
#define MAX_STATES   8
#define MAX_NODES    8

#define TOGPU 0
#define FROMGPU 1

#define MAXPRODUCT 1
#define SUMPRODUCT 0

typedef struct
{
  int __align__(8) size; /*number of states*/
  int __align__(8) npar; /*number of parents*/
  int __align__(8) parents[MAX_PARENTS]; 
  int __align__(8) nchild;
  int __align__(8) children[MAX_CHILDREN];
  int __align__(8) observed;
  float __align__(8) *M; 
  int   __align__(8) Msize;
  float __align__(8) lambda[MAX_STATES];
  float __align__(8) pi[MAX_STATES];
  float __align__(8) bel[MAX_STATES];
}node;


typedef struct __tagmessage
{
  float __align__(8) lambda [MAX_NODES][MAX_PARENTS][MAX_STATES];
  float __align__(8) pi     [MAX_NODES][MAX_CHILDREN][MAX_STATES];
}message;


typedef struct __tagnetwork
{
  int __align__(8) num;
  int __align__(8) sizes[MAX_NODES];
  node __align__(8) var[MAX_NODES];
  int __align__(8) dag[MAX_NODES][MAX_PARENTS];
}network;

typedef struct __tagstate
{
  float __align__(8) value[MAX_NODES][MAX_STATES];
  int   __align__(8) sizes[MAX_NODES];
  int   __align__(8) num;
}state;

typedef struct
{
  int  __align__(8) cidx[MAX_CHILDREN];
  int  __align__(8) pidx[MAX_PARENTS];
  int  __align__(8) psizes[MAX_PARENTS];
  int  __align__(8) uidx[MAX_PARENTS];
  int  __align__(8) msize;
  float __align__(8) lambda[MAX_STATES];
  float __align__(8) pi[MAX_STATES]; 
  float __align__(8) bel[MAX_STATES];
  float __align__(8) term[MAX_NODES][MAX_STATES][MAX_STATES];
  float __align__(8) py[MAX_STATES];
}frame;

#endif 
