function tests=test_b1_forward_bi_p2c_gain
tests=functiontests(localfunctions);
end

function testEnergyOddResponse(test)
B=struct('E',(100:100:800).','q',[-.002 -.001 .001 .002]);
B.Y=zeros(8,4,3);
[X,lim]=b1_forward_bi_p2c_gain(B,'gainE');
verifySize(test,X,[8 4 3 6]); verifyNumElements(test,lim,6);
verifyEqual(test,X(:,1,:,:),-X(:,4,:,:));
verifyEqual(test,X(:,2,:,:),-X(:,3,:,:));
verifyEqual(test,X(1,:,:,:),X(2,:,:,:));
verifyEqual(test,X(7,:,:,4:6),zeros(1,4,3,3));
[X,lim]=b1_forward_bi_p2c_gain(B,'none');
verifyEmpty(test,lim); verifySize(test,X,[8 4 3 0]);
end

function testShiftDerivativeSign(test)
B=struct('E',(100:100:600).','q',[-.002 -.001 .001 .002]);
B.Y=zeros(6,4,3);
B.archive=repmat(exp(9-600*abs(B.q)),6,1,3);
B.intr=qe_lowq_joint_forward('intrinsic');
[X,~]=b1_forward_bi_p2c_gain(B,'shift');
verifyEqual(test,X(:,3,:,1),.6*ones(6,1,3),'AbsTol',1e-12);
verifyEqual(test,X(:,3,:,1),-X(:,2,:,1),'AbsTol',1e-12);
end
