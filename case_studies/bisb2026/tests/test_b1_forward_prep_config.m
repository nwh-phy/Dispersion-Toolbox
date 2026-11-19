function tests = test_b1_forward_prep_config
tests = functiontests(localfunctions);
end

function testBiAcquisitionsUseSixContiguousBlocks(t)
ids = {'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'};
for k = 1:numel(ids)
    cfg = b1_forward_prep_config(ids{k},tempname);
    d = bisb_datasets(ids{k});
    verifyEqual(t,cfg.n_blocks,6);
    verifyEqual(t,cfg.frames_per_block,d.n_frames/6);
    verifyGreaterThanOrEqual(t,cfg.qmax,b1_p2a_inputs_config(cfg.paths).q_max);
end
end
