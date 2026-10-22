function tests = test_browser_session_io
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testStreamingHashMatchesStandardBinaryDigest(testCase)
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
path = fullfile(folder, 'binary.dat');
fid = fopen(path, 'wb'); fwrite(fid, uint8(0:255)); fclose(fid);
verifyEqual(testCase, qe_file_sha256(path), ...
    '40aff2e9d2d8922e47afd4648e6967497158785fbd1da870e7110266bf944880');
end

function testSourceIdentityDetectsSameSizeChanges(testCase)
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
path = fullfile(folder, 'eq3D.mat');
fid = fopen(path, 'wb'); fwrite(fid, uint8('abc')); fclose(fid);
sources = qe_browser_source_identity(struct('source_path', path));
verifyEmpty(testCase, qe_check_browser_sources(sources));
fid = fopen(path, 'wb'); fwrite(fid, uint8('abd')); fclose(fid);
verifyError(testCase, @() qe_check_browser_sources(sources), 'qe_check_browser_sources:SourceChanged');
delete(path);
verifyEqual(testCase, qe_check_browser_sources(sources), {path});
end

function testDataFingerprintBindsCountsAndCalibration(testCase)
qe = make_qe_struct(reshape(1:60, 10, 6), (1:10)', 0.0005);
baseline = qe_data_fingerprint(qe);
changed = qe; changed.intensity(3) = changed.intensity(3) + 1;
verifyNotEqual(testCase, qe_data_fingerprint(changed), baseline);
changed = qe; changed.energy_meV = changed.energy_meV + 0.01;
verifyNotEqual(testCase, qe_data_fingerprint(changed), baseline);
changed = qe; changed.dq_Ainv = 0.00025;
verifyNotEqual(testCase, qe_data_fingerprint(changed), baseline);
end
