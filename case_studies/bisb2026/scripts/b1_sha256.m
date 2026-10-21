function hex = b1_sha256(x)
% SHA-256 (lower-case hex) of a char/string (UTF-8 bytes), a uint8 array, a
% struct (hashed through jsonencode), or the bytes of a file given as
% {'file', path}. Used as configuration / code fingerprint in run folders.
if iscell(x) && numel(x)==2 && strcmp(x{1},'file')
 fid=fopen(x{2},'r'); assert(fid>0,'b1_sha256:File','cannot open %s',x{2}); b=fread(fid,'*uint8'); fclose(fid);
elseif isstruct(x)
 b=unicode2native(jsonencode(x),'UTF-8');
elseif ischar(x) || isstring(x)
 b=unicode2native(char(x),'UTF-8');
else
 b=uint8(x);
end
md=java.security.MessageDigest.getInstance('SHA-256');
h=typecast(md.digest(b(:)),'uint8');
hex=lower(reshape(dec2hex(h,2).',1,[]));
end
