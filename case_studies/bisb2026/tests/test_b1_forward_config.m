function tests = test_b1_forward_config
tests = functiontests(localfunctions);
end

function testLegacyLocationsAndFingerprints(t)
p = b1_forward_config();
verifyEqual(t,p.session,'590_PL2_10w');
verifyTrue(t,endsWith(p.input_dir,fullfile('p2a','inputs_590_PL2_10w_66385d96')));
verifyEqual(t,p.response_dir,fullfile(p.input_dir,'response_42ffc295'));
verifyEqual(t,b1_sha256(b1_p2b_config()), ...
    '691c36adfee10ee270c8b15c7536be4e7b8c5730362478f82233132fd8299a58');
[c,o] = b1_p2b_slot_config('h66');
verifyEqual(t,b1_sha256(c),'3fb982f21aaf84180895f19c11f3f452758566c18f19c5e86ab2c25c4f273880');
verifyEqual(t,[o.h o.n_sy],[0.0066 53]);
end

function testConfigured20wWithAdditionalNode(t)
p = b1_forward_config('no_PL2_20w_2film', tempname);
n = b1_p2b_config().nodes;
n.E0 = sort([n.E0 0.0015]);
cfg = b1_p2b_config(p,n);
verifyEqual(t,cfg.nodes.E0,n.E0);
verifySize(t,cfg.starts(1).E0,[2 6]);
verifySize(t,cfg.starts(1).G,[2 3]);
verifyEqual(t,numel(cfg.exit_q),40);
verifyEqual(t,min(abs(cfg.exit_q)),bisb_datasets(p.session).dq);
verifyEqual(t,p.fit_prefix,'fit_no_PL2_20w_2film_');
[c,o] = b1_p2b_slot_config('hscan',0.0033,cfg);
verifyEqual(t,c.nodes,n);
verifyEqual(t,o.h,0.0033);
verifyEqual(t,o.tail_edges(1),0.006625,'AbsTol',1e-12);
end
