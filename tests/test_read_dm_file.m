function tests=test_read_dm_file
tests=functiontests(localfunctions);
end
function testFourDimensionalFirstImage(t)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'third_party','read_dm_file'));
p=[tempname '.dm4']; clean=onCleanup(@()delete(p));
expected=reshape(single(1:24),[2 3 2 2]); write_fixture(p,expected,1);
[m,s,u,meta]=ReadDMFile(p,[p '.log'],1); delete([p '.log']);
verifyEqual(t,double(m),double(expected)); verifyEqual(t,double(s),0.004,'AbsTol',1e-9); verifyEqual(t,u,'eV'); verifyEqual(t,meta{8},single(100));
end
function testOriginalSecondImageDefault(t)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'third_party','read_dm_file'));
p=[tempname '.dm4']; clean=onCleanup(@()delete(p));
expected=reshape(single(1:12),[2 3 2]); write_fixture(p,expected,2);
[m,s,u]=ReadDMFile(p); verifyEqual(t,double(m),double(expected)); verifyEqual(t,double(s),0.004,'AbsTol',1e-9); verifyEqual(t,u,'eV');
end
function write_fixture(path,A,index)
dims={}; sz=size(A);
for k=1:numel(sz), dims{end+1}=leaf('',[3],typecast(int32(sz(k)),'uint8')); end
cal=group('Calibrations',{group('Dimension',{group('',{leaf('Scale',6,typecast(single(.004),'uint8')),leaf('Units',[20 4 2],typecast(uint16('eV'),'uint8')),leaf('Origin',6,typecast(single(100),'uint8'))})})});
im=group('',{group('ImageData',{group('Dimensions',dims),leaf('Data',[20 6 numel(A)],typecast(A(:)','uint8')),cal})});
images=repmat({group('',{})},1,index); images{index}=im;
body=uint8([1 1 be64(1) group('ImageList',images)]);
f=fopen(path,'w'); fwrite(f,[be32(4) be64(numel(body)) be32(1) body],'uint8'); fclose(f);
end
function b=group(name,children)
payload=uint8([1 1 be64(numel(children)) children{:}]);
b=uint8([20 be16(numel(name)) uint8(name) be64(numel(payload)) payload]);
end
function b=leaf(name,info,data)
v=[]; for x=info, v=[v be64(x)]; end
payload=uint8([uint8('%%%%') be64(numel(info)) v data]);
b=uint8([21 be16(numel(name)) uint8(name) be64(numel(payload)) payload]);
end
function b=be16(x), b=typecast(swapbytes(uint16(x)),'uint8'); end
function b=be32(x), b=typecast(swapbytes(uint32(x)),'uint8'); end
function b=be64(x), b=typecast(swapbytes(uint64(x)),'uint8'); end
