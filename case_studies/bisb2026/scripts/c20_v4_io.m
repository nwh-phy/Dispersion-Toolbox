function result = c20_v4_io(action,path,value)
result=[];
switch action
 case 'hash'
  md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(path,'rb'); assert(fid>=0);
  c=onCleanup(@()fclose(fid));
  while ~feof(fid), x=fread(fid,8*1024*1024,'*uint8'); md.update(typecast(x,'int8')); end
  result=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
 case 'digest'
  md=java.security.MessageDigest.getInstance('SHA-256'); md.update(unicode2native(path,'UTF-8'));
  result=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
 case {'json','text'}
  if strcmp(action,'json'), value=jsonencode(value,PrettyPrint=true); end
  fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',value);
 case 'inventory'
  files=dir(fullfile(path,'**','*')); files=files(~[files.isdir]); rows=cell(numel(files),3);
  for k=1:numel(files)
   p=fullfile(files(k).folder,files(k).name);
   rows(k,:)={string(extractAfter(p,[char(path) filesep])),files(k).bytes,string(c20_v4_io('hash',p))};
  end
  result=cell2table(rows,'VariableNames',{'path','bytes','sha256'});
end
end
