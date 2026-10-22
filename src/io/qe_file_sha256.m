function hash = qe_file_sha256(path)
%QE_FILE_SHA256 Stream a file into SHA-256 without loading it into memory.
fid = fopen(path, 'rb');
if fid < 0
    error('qe_file_sha256:OpenFailed', 'Cannot read source file: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');
fseek(fid, 0, 'eof');
remaining = ftell(fid);
fseek(fid, 0, 'bof');
while remaining > 0
    count = min(remaining, 4 * 1024 * 1024);
    block = fread(fid, count, '*uint8');
    if numel(block) ~= count
        error('qe_file_sha256:ReadFailed', 'Could not finish reading: %s', path);
    end
    md.update(typecast(block, 'int8'));
    remaining = remaining - count;
end
hash = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
