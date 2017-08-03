function [base,ext]=fileparts2(string)
[path,name,ext]=fileparts(deblank(string));
if strcmp(ext,'.gz')
   [path2,name,ext]=fileparts(name);
end   
if isempty(path)
   base=name;
else
   base=[path '/' name];
end
return

