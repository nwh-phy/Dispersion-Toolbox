function tests = test_b1_forward_bi_p2b_material
tests = functiontests(localfunctions);
end

function testPreservesDonorInterpolationOnCommonSupport(t)
c = struct('nodes',struct('E0',[.00025 .001 .005], ...
    'G',[.00025 .002 .005],'A',[.00025 .002 .005]),'variants',struct('ext','const'));
e = [500 700 800;1100 1200 1300];
g = [300 500 800;500 1000 1800];
a = [0.3 1 2;1 0.6 1];
f = struct('E0',e,'G',g,'a',reshape(a.',[],1), ...
    'u',[reshape(log(e.^2).',1,[]),reshape(log(g).',1,[])]);
S = struct('cfg',c,'res',struct('fit',f));
a = b1_forward_bi_p2b_material(S);
b = b1_forward_bi_p2b_material(S,.0005);
q = [.0005 .0006 .0008 .001 .002 .003 .005];
e = 250:4:1800;
L = qe_lowq_joint_forward('C',q,e,a.theta,a.intrinsic);
R = qe_lowq_joint_forward('C',q,e,b.theta,b.intrinsic);
L = sum(L.*reshape(a.a,1,1,[]),3);
R = sum(R.*reshape(b.a,1,1,[]),3);
verifyEqual(t,L,R,'RelTol',5e-13);
verifyEqual(t,b.intrinsic.Qn_E0(1),.0005);
verifyEqual(t,b.intrinsic.Qn_G(1),.0005);
verifyEqual(t,b.intrinsic.Qn_A(1),.0005);
end
