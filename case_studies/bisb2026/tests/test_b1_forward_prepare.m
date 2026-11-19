function tests = test_b1_forward_prepare
tests = functiontests(localfunctions);
end

function testRejectExistingP1bBeforeReadingRaw(t)
folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
cfg = b1_forward_prep_config('Bi1107_10w',folder.Folder);
mkdir(cfg.paths.diag_dir);
verifyError(t,@() b1_forward_prepare(cfg),'b1_forward_prepare:OutputExists');
verifyFalse(t,isfile(cfg.paths.stage1_file));
end
