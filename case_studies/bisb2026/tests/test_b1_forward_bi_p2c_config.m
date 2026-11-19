function tests=test_b1_forward_bi_p2c_config
tests=functiontests(localfunctions);
end

function testWindowsAndRegisteredGeometry(test)
ids={'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'};
lo=[200 200 150 250];
ctrl=[180 180 180 300];
sy=[109 57 109 57];
for i=1:numel(ids)
    c=b1_forward_bi_p2c_config(ids{i},'unit_test');
    verifyEqual(test,c.window,[lo(i) 1800]);
    verifyEqual(test,c.control,[ctrl(i) 1800]);
    verifyEqual(test,c.nsy,sy(i));
    verifyEqual(test,c.o.h,bisb_datasets(ids{i}).h);
end
c=b1_forward_bi_p2c_config('Bi1107_20w','unit_test');
verifyEqual(test,c.qmax,12*bisb_datasets(c.id).dq);
end
