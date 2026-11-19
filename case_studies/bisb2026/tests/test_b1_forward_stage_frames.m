function tests = test_b1_forward_stage_frames
tests = functiontests(localfunctions);
end

function testCenterRestoredInNativeOrder(t)
folder = t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
file = fullfile(folder.Folder,'channels.mat');
C.q = [-2 -1 1 2] * 0.00025;
frames = reshape(single(1:24),3,2,4);
center_frames = single([12 20;91 7;12 13]);
save(file,'frames','center_frames','-v7.3');
[A,q] = b1_forward_stage_frames(file,C,0.00025);
verifyEqual(t,q,[-1 0 1]*0.00025);
verifyEqual(t,A,permute(double(cat(3,frames(:,:,2),center_frames,frames(:,:,3))),[1 3 2]));
end
