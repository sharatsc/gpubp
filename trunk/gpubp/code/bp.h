#include "common.h"
#include <iostream>
#include <fstream>
#include "cuda.h"


using namespace std;
extern message h_src,h_dst; /*host messages*/
extern frame h_frame[MAX_NODES]; /*host frames*/

/*output function*/
ostream & operator<<(ostream & out,const network& net);
istream & operator>>(istream& in,network& net);
ostream& operator<<(ostream& out,const state& st);
istream& operator>>(istream& in,state& st);


/*network functions*/
void free_network_resources(network* pnet);
void getidx(const int *psizes,int* uidx,int npar,int m);
void init_messages(message* pm);



