function tests = test_qe_align_sequence_frames
tests = functiontests(localfunctions);
end

function testDetectorAndReferenceEnergyZero(t)
folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
source = fullfile(folder.Folder,'counts.npy');
raw = ones(3,9,51,'uint32');
for frame = 1:3
    raw(frame,2:8,24+frame) = uint32([10 25 45 70 50 30 15]);
end
raw(:,9,28) = 1000;                  % global detector sum peaks at 28
raw(2,1,32) = intmax('uint32');      % ignore the entire first column in the summary
write_npy_fixture(source,raw);
metadata.spatial_calibrations = struct('scale',{1,1,0.004},'offset',{0,0,-0.1},'units',{'','','eV'});
writelines(jsonencode(metadata),fullfile(folder.Folder,'counts.json'));
[A,R] = qe_align_sequence_frames(source,5,0.0005,q_select=0.0015,energy_zero="detector_peak");
B = qe_align_sequence_frames(source,5,0.0005,q_select=0.0015);
[B2,R2] = qe_align_sequence_frames(source,5,0.0005,q_select=0.0015);
verifyEqual(t,B,B2);
verifyEqual(t,R,R2);
verifyEqual(t,R.valid_channels,2:9);
verifyEqual(t,R.peak_pixel,28);
verifyEqual(t,R.counts_by_energy_frame,squeeze(sum(double(raw(:,2:9,:)),2)).');
verifyEqual(t,A.alignment.measured_offset_pixels,[-1 0 1]);
verifyEqual(t,A.alignment.support,(2:50).');
verifyEqual(t,A.E,B.E-8,'AbsTol',1e-10);
verifyEqual(t,A.F,B.F);
end

function write_npy_fixture(source,raw)
shape = size(raw);
header = sprintf('{''descr'': ''<u4'', ''fortran_order'': True, ''shape'': (%d, %d, %d), }',shape);
padding = mod(-10-numel(header)-1,16);
header = [header repmat(' ',1,padding) newline];
fid = fopen(source,'w','ieee-le');
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid,[147 uint8('NUMPY') 1 0],'uint8');
fwrite(fid,numel(header),'uint16');
fwrite(fid,header,'char');
fwrite(fid,raw,'uint32');
end
