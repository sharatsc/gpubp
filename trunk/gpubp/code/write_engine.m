%-----------------------------------------------------
%format for the network
%num nodes
%
%sharat@mit.edu
%
%
function write_engine(engine,destfile)
   fid   =fopen(destfile,'w');
   engine=struct(engine);
   bnet  =struct(engine.inf_engine);
   bnet  =struct(bnet.bnet);
   dag   =bnet.dag;bnet.dag=[];
   num_nodes=size(dag,1);
   fprintf(fid,'%d\n',num_nodes);%num nodes
   for n=1:num_nodes
	 parents=find(dag(:,n));
	 children=find(dag(n,:));
	 cpt     =[];
	 cpd     =struct(bnet.CPD{n});
	 if(~isempty(cpd))
	   cpt     =cpd.CPT;
	 end;
	 fprintf(fid,'#-----------------node(%d)\n',n);
	 fprintf(fid,'#size\n');
	 fprintf(fid,'%d\n',bnet.node_sizes(n));
	 fprintf(fid,'#parents\n',n);
	 fprintf(fid,'%d ',[length(parents)]);
	 for p=1:length(parents)
	   fprintf(fid,'%d ',parents(p));
	 end;
     fprintf(fid,'\n');
	 fprintf(fid,'#children\n',n);
	 fprintf(fid,'%d ',[length(children)]);
	 for c=1:length(children)
	   fprintf(fid,'%d ',children(c));
	 end;
     fprintf(fid,'\n');
	 fprintf(fid,'#M\n');
	 fprintf(fid,'%d ',length(cpt(:)));
	 fprintf(fid,'%.6f ',cpt(:));
	 fprintf(fid,'\n');
   end;
   fclose(fid);
%end function
