%--------------------------------------------------------------------------------------
%
%sharat@mit.edu
%
function matlab_baseline(nodes)
   HOME=''; %set to root of BNT source folder
   if(strcmp(HOME,''))
	 HOME=pwd()
   end;
   %addpath(genpath(HOME));
   %----------------------------------------
   %set up network
   if(1)
	 dag             = create_simple_dag;
	 sizes           = 2*[1 1 1];
	 parents{1}      = [];
	 parents{2}      = 1;
	 parents{3}      = 1;
   elseif(0)
	 levels          = ceil(log2(nodes));
     [dag,parents]   = create_tree_dag(levels);
	 sizes           = 2*ones(size(dag,1),1);%+floor(5*rand(1,size(dag,1)));
   elseif(1)
     [dag,parents]   = create_polytree_dag(nodes);
	 sizes           = 2+round(2*rand(1,size(dag,1)));
   else
     [dag,parents,sizes]   = create_inf_dag(nodes);
   end;
   %-----------
   %sizes
   bnet            = mk_bnet(dag,sizes);
   dsize           = size(dag,1);
   %-----------
   %set up CPT
   for i=1:dsize
	  if(isempty(parents{i}))
	    bnet.CPD{i}=tabular_CPD(bnet,i,'CPT','unif');
	  else;
		sz         =sizes(i);
	    psz        =prod(sizes(parents{i}));
        cpt        =rand(psz,sz);cpt=cpt./repmat(sum(cpt,2),1,sz);
	    bnet.CPD{i}=tabular_CPD(bnet,i,'CPT',cpt);
	 end;
   end;
   engine=jtree_inf_engine(bnet);%,'protocol','parallel','filename','/dev/null');
   %----------------------------------------
   %generate input
   evidence=cell(dsize,1);
   sevidnce=cell(dsize,1);
   output  =cell(dsize,1);
   for i=1:dsize
	 if((parents{i}))%child node
	   lambda      =rand(sizes(i),1);
	   sevidence{i}=lambda/sum(lambda);
	 end;
   end;
   write_engine(engine,'net.txt');
   write_state(sevidence,'net-input.txt');
   
   tic;
   engine=enter_evidence(engine,evidence,'soft',sevidence);
   for i=1:dsize
	 marg=marginal_nodes(engine,i);
	 output{i}=marg.T(:);
   end;
   time=toc;
   
   
   write_state(output,'matlab-output.txt');
   fprintf('MATLAB time elapsed for %d nodes is:%.6f sec\n',nodes,time);
%end function

function [dag,parents]=create_simple_dag
   dag=zeros(3);
   sizes=[2 2 2];
   dag(1,2)=1;
   dag(1,3)=1;
   
function [dag,parents,sizes]=create_inf_dag(nodes)
   NLOC=16;
   NOBJ=5;
   
   N=NLOC^2;
   L=1;O=2;
   
   F_start=2;
   
   C_start=F_start+nodes;
   dsize=C_start+nodes;
   dag=zeros(dsize);
   sizes=[N NOBJ 2*ones(1,nodes) (N+1)*ones(1,nodes)];
   
   for n=1:nodes
	 parents{F_start+n}=[O];
	 parents{C_start+n}=[L F_start+n];
	 dag(O,F_start+n)        =1;
	 dag(F_start+n,C_start+n)=1;
	 dag(L,C_start+n)        =1;
   end;
   
function [dag,parents]=create_polytree_dag(dsize)
   idx=randperm(dsize);n0=ceil(dsize/3);n1=floor(dsize/3);n2=dsize-n0-n1;
   l0 =idx(1:n0);idx(1:n0)=[];
   l1 =idx(1:n1);idx(1:n1)=[];
   l2 =idx(1:n2);idx(1:n2)=[];
   idx
   dag=zeros(dsize);
   parents=cell(1,dsize);
   for i=1:length(l1)
	 parents{l1(i)}=l0;
	 dag(l0,l1(i))=1;
   end;
   for i=1:length(l2)
	 parents{l2(i)}=l1;
	 dag(l1,l2(i))=1;
   end;

%------------------------------
%creates a tree with levels
function [dag,parents]=create_tree_dag(levels)
    curparents ={1};
	nextparents={1};
	parents{1} =[];
	idx        = 2;
	for l=2:levels
	  curparents=nextparents;
	  nextparents={};nd=1;
	  for p=1:length(curparents)
		parents{idx}=curparents{p};nextparents{nd}=idx;idx=idx+1;nd=nd+1;
		parents{idx}=curparents{p};nextparents{nd}=idx;idx=idx+1;nd=nd+1;
      end;
    end;
	dag=zeros(length(parents));
	for n=1:length(parents)
	  for p=1:length(parents{n})
		dag(parents{n}(p),n)=1;
	  end;
	end;
%end function
