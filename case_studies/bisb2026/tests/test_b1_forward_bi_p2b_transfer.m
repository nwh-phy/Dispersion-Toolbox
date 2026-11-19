function tests = test_b1_forward_bi_p2b_transfer
tests = functiontests(localfunctions);
end

function testRecoversScaleAndFreezesItBeforePrediction(t)
[P,F,c] = fixture();
X = b1_forward_bi_p2b_transfer(P,F,c);
verifyEqual(t,X.scale,2,'AbsTol',1e-14);
verifyEqual(t,X.table.cross_chi2(1),0,'AbsTol',1e-14);
verifyEqual(t,X.table.cross_chi2(2),X.table.own_chi2(2),'AbsTol',1e-14);
verifyEqual(t,X.table.cross_dof,[3;4]);
P.Y(:,:,2) = 999;
Z = b1_forward_bi_p2b_transfer(P,F,c);
verifyEqual(t,Z.scale,X.scale);
verifyEqual(t,Z.mu,X.mu);
verifyGreaterThan(t,Z.table.cross_chi2(2),X.table.cross_chi2(2));
end

function testNonnegativeScaleAndSplitValidation(t)
[P,F,c] = fixture();
P.Y(:,:,1) = 0;
X = b1_forward_bi_p2b_transfer(P,F,c);
verifyEqual(t,X.scale,0);
F.mask(:,:,2) = true;
verifyError(t,@() b1_forward_bi_p2b_transfer(P,F,c),'b1_forward_bi_p2b_transfer:Split');
end

function [P,F,c] = fixture()
c = struct('holdout',[150 250],'lobo_block',2,'tail',struct('E_ref',250,'p',4));
P = struct('E',[250;350],'mu',reshape(1:8,2,2,2),'S2',reshape(3:10,2,2,2),'mask',true(2,2,2));
F = struct('mask',P.mask,'b',[1;2;3;4],'tailOn',true,'dof',2);
F.mask(:,:,2) = false;
tail = (P.E/250).^(-4).*reshape(F.b,1,2,2);
F.mu = 2*P.mu+tail;
P.Y = F.mu;
P.Y(:,:,2) = F.mu(:,:,2)+3;
end
