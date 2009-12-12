%--------------------------------------------
%
%sharat@mit.edu
%
function state=read_state(filename)
   fid=fopen(filename,'r');
   nstates=fscanf(fid,'%d',1)
   state=cell(nstates,1);
   for i=1:nstates
	 state{i}=[];
	 n=fscanf(fid,'%d',1)
	 size=fscanf(fid,'%d',1)
	 state{n}=fscanf(fid,'%f',size)
   end;
   fclose(fid);
%end function
