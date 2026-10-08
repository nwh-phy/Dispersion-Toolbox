function result = c20_packet_readback(archive,report,verify)
% Independently extract a bounded packet, persist verification, then clean up.
assert(~isfile(report),'c20:ReportExists','Readback report exists; refusing overwrite');
zip_hash=c20_v4_io('hash',archive);
base=char(java.io.File(tempdir).getCanonicalPath());
attrs=javaArray('java.nio.file.attribute.FileAttribute',0);
decoded=char(java.nio.file.Files.createTempDirectory(java.io.File(base).toPath(),'c20-readback-',attrs).toString());
try
    z=java.util.zip.ZipFile(archive); close_zip=onCleanup(@()z.close());
    entries=z.entries(); names={}; bytes=0;
    while entries.hasMoreElements()
        e=entries.nextElement(); name=char(e.getName());
        safeChild(decoded,name);
        assert(~ismember(name,names),'c20:DuplicateEntry','Duplicate ZIP entry');
        names{end+1}=name; %#ok<AGROW>
        assert(e.getSize()>=0,'c20:UnknownSize','ZIP entry size unknown');
        bytes=bytes+double(e.getSize());
        assert(bytes<=1024^3 && numel(names)<=10000,'c20:PacketLimit','Packet exceeds readback diagnostic limit');
    end
    clear close_zip
    unzip(archive,decoded);
    manifest=readtable(fullfile(decoded,'FILE_MANIFEST.csv'),TextType='string');
    files=c20_v4_io('inventory',decoded);
    paths=replace(manifest.path,char(92),'/');
    actual=replace(files.path,char(92),'/'); actual(actual=="FILE_MANIFEST.csv")=[];
    assert(numel(unique(paths))==height(manifest) && isequal(sort(paths),sort(actual)), ...
        'c20:FileSet','Manifest does not match complete payload set');
    for k=1:height(manifest)
        file=safeChild(decoded,char(manifest.path(k))); info=dir(file);
        assert(isscalar(info) && info.bytes==manifest.bytes(k) && ...
            strcmpi(c20_v4_io('hash',file),manifest.sha256(k)), ...
            'c20:HashMismatch','ZIP payload hash or size mismatch');
    end
    result=verify(decoded);
    assert(strcmpi(zip_hash,c20_v4_io('hash',archive)),'c20:ZipChanged','ZIP changed during readback');
    summary=struct('archive',archive,'zip_sha256',zip_hash,'file_count',height(files), ...
        'payload_hashes',height(manifest),'passed',true,'temporary_directory',decoded, ...
        'cleanup_status','pending','verification',result);
    c20_v4_io('json',report,summary);
catch err
    warning('c20:ReadbackFailed','Readback failed; diagnostic directory retained: %s',decoded);
    rethrow(err);
end
try
    safeChild(base,char(java.io.File(decoded).getName()));
    assert(startsWith(char(java.io.File(decoded).getName()),'c20-readback-'));
    assert(~java.nio.file.Files.isSymbolicLink(java.io.File(decoded).toPath()));
    rmdir(decoded,'s'); summary.cleanup_status='removed';
catch err
    summary.cleanup_status='failed'; summary.cleanup_error=err.message;
    warning('c20:CleanupFailed','Readback cleanup failed; directory retained: %s (%s)',decoded,err.message);
end
c20_v4_io('json',report,summary);
end

function path=safeChild(base,name)
name=strrep(name,char(92),'/');
if startsWith(name,'/') || ~isempty(regexp(name,'(^|/)\.\.(/|$)|:','once'))
    error('c20:PathEscape','Unsafe packet path: %s',name);
end
base=char(java.io.File(base).getCanonicalPath());
path=char(java.io.File(base,name).getCanonicalPath());
assert(startsWith(lower(path),[lower(base) filesep]),'c20:PathEscape','Path outside owned directory');
end
