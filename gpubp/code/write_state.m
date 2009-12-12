%-----------------------------------------------
%
%sharat@mit.edu
%-------------------------------------------------
function write_state(states,filename)
fid=fopen(filename,'w');
num=length(states)
fprintf(fid,'%d\n',num);
for i=1:length(states)
  fprintf(fid,'%d %d ',i,length(states{i}));
  fprintf(fid,'%.6f ',states{i});
  fprintf(fid,'\n');
end;  
fclose(fid);
%end function
