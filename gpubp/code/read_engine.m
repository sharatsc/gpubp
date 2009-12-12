%-----------------------------------------------------
%format for the network
%num nodes
%
%sharat@mit.edu
%
%
function engine=read_engine(srcfile)
   fid   =fopen(srcfile,'r');
   %read num_nodes
   str   =read_string(fid);
   num_nodes=sscanf(str,'%d',[1]);
   dag   =zeros(num_nodes);
   sizes =zeros(num_nodes,1);
   CPD   ={};
   for n=1:num_nodes
	 fprintf('reading node:%d\n',n);
	 %size
	 str=read_string(fid);
	 sizes(n)=sscanf(str,'%d');
	 
	 %parents
	 str=read_string(fid);
	 parents=sscanf(str,'%d');
	 npar   =parents(1);
     parents=parents(2:end);
	 if(parents)
	   dag(parents,n)=1;
	 end;

	 %children
	 str=read_string(fid);
	 children=sscanf(str,'%d');
	 nchild   =children(1);
     children=children(2:end);
	 if(children)
	   dag(n,children)=1;
	 end;
	 %cpt
	 str=read_string(fid);
	 cpt=sscanf(str,'%f');
	 CPD{n}=cpt(2:cpt(1)+1);
   end;
   %build whole model
   bnet=mk_bnet(dag,sizes);
   for n=1:num_nodes
	 fprintf('building:%d\n',n);
	 bnet.CPD{n}=tabular_CPD(bnet,n,'CPT',CPD{n}(:));
   end;
   engine=jtree_inf_engine(bnet);
   fclose(fid);
%end function

function str=read_string(fid)
  str='';
  while(1)
	str=fgetl(fid);
	if(isnumeric(str)) break;end;
	if(str(1)=='#') continue;end;
	break;
  end;
%end function
