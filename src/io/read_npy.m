function data = read_npy(filename)
%READ_NPY  Read a NumPy .npy file into a MATLAB array.
%
%   data = read_npy('file.npy')
%
%   Supports N-D arrays of standard numeric types (uint8/16/32/64,
%   int8/16/32/64, single, double, logical).  Both C-order and
%   Fortran-order layouts are handled.
%
%   Based on npy-matlab (https://github.com/kwikteam/npy-matlab).
%   Bundled here for self-contained use without external dependencies.
%
%   See also: load_npy_session, load_qe_dataset

arguments
    filename {mustBeFile}
end

%% Parse header
[arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength] = ...
    local_parse_header(filename);

%% Read data
if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try
    fread(fid, totalHeaderLength, 'uint8');   % skip header
    data = fread(fid, prod(arrayShape), [dataType '=>' dataType]);

    if length(arrayShape) > 1 && ~fortranOrder
        data = reshape(data, arrayShape(end:-1:1));
        data = permute(data, length(arrayShape):-1:1);
    elseif length(arrayShape) > 1
        data = reshape(data, arrayShape);
    end

    fclose(fid);
catch me
    fclose(fid);
    rethrow(me);
end
end


function [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength] = ...
        local_parse_header(filename)
%LOCAL_PARSE_HEADER  Parse the binary header of a .npy file.

dtypesMatlab = {'uint8','uint16','uint32','uint64', ...
                'int8','int16','int32','int64', ...
                'single','double','logical'};
dtypesNPY    = {'u1','u2','u4','u8', ...
                'i1','i2','i4','i8', ...
                'f4','f8','b1'};

fid = fopen(filename);
if fid == -1
    if ~isempty(dir(filename))
        error('read_npy:PermissionDenied', 'Permission denied: %s', filename);
    else
        error('read_npy:FileNotFound', 'File not found: %s', filename);
    end
end

try
    % Magic string: \x93NUMPY
    magicString = fread(fid, [1 6], 'uint8=>uint8');
    if ~all(magicString == [147, 78, 85, 77, 80, 89])
        error('read_npy:NotNPY', 'Not a valid .npy file: %s', filename);
    end

    majorVersion = fread(fid, [1 1], 'uint8=>uint8');  %#ok<NASGU>
    minorVersion = fread(fid, [1 1], 'uint8=>uint8');  %#ok<NASGU>
    headerLength = fread(fid, [1 1], 'uint16=>uint16');
    totalHeaderLength = 10 + headerLength;

    arrayFormat = fread(fid, [1 headerLength], 'char=>char');

    % Data type
    r = regexp(arrayFormat, '''descr''\s*:\s*''(.*?)''', 'tokens');
    dtNPY = r{1}{1};
    littleEndian = ~strcmp(dtNPY(1), '>');
    dataType = dtypesMatlab{strcmp(dtNPY(2:3), dtypesNPY)};

    % Fortran order
    r = regexp(arrayFormat, '''fortran_order''\s*:\s*(\w+)', 'tokens');
    fortranOrder = strcmp(r{1}{1}, 'True');

    % Shape
    r = regexp(arrayFormat, '''shape''\s*:\s*\((.*?)\)', 'tokens');
    shapeStr = r{1}{1};
    arrayShape = str2num(shapeStr(shapeStr ~= 'L'));  %#ok<ST2NM>

    fclose(fid);
catch me
    fclose(fid);
    rethrow(me);
end
end
