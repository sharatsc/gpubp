%---------------------------------------------------------------------------
%
%sharat@mit.edu
%---------------------------------------------------------------------------
HOME=''; %set to root of BNT source folder
if(strcmp(HOME,''))
  HOME=pwd()
end;
addpath(genpath(fullfile(HOME,'code')));
addpath(genpath(fullfile(HOME,'BNT')));
