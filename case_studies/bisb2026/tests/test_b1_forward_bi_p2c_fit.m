function tests=test_b1_forward_bi_p2c_fit
tests=functiontests(localfunctions);
end

function testLegacyModePredictionAndMask(test)
q=[-.002 -.001 .001 .002];
Q=[0;.001;.0015;.002;.003];
W=zeros(4,1,5);
for i=1:4
    j=2+(abs(q(i))>.001)*2;
    W(i,1,j)=1e-5;
end
ops=repmat({struct('W',W,'Qm',Q)},1,6);
C=qe_lowq_joint_forward('build',ops,repmat({[0 1 0]},1,6), ...
    struct('omega',-200:4:1800,'eps',-4:4:4,'variant','slice'));
nodes=struct('E0',[.001 .002 .003],'G',[.001 .003],'A',[.001 .003]);
I=qe_lowq_joint_forward('intrinsic','Qn_E0',nodes.E0,'Qn_G',nodes.G,'Qn_A',nodes.A);
th=struct('E0sq',[450 560 640;1100 1180 1300].^2,'G',[210 350;390 600]);
a=qe_lowq_joint_forward('eval',C,th,I);
y=reshape(sum(a.*reshape([.4 .6 .3 .8]*1e8,1,[],1,1),2),numel(C.E_out),4,6);
tail=max(C.E_out(:)/250,150/250).^-4; tail(C.E_out<150)=0;
y=y+tail.*reshape(3*(1:4),1,4,1);
cfg=struct('nodes',nodes,'bounds_E0',[150 950;700 2500],'bounds_G',[20 2400], ...
    'tail',struct('E_ref',250,'p',4),'holdout',[150 250],'lambda_curv',0);
B=struct('E',C.E_out(:),'Y',y,'S2',max(y,1),'q',q,'cfg',cfg,'intr',I, ...
    'cache',C,'mask',true(size(y)),'exposure',ones(6,1),'stage','odd');
seed=struct('u',[reshape(log(th.E0sq).',1,[]) reshape(log(th.G).',1,[])]);
P=struct('gain','none','window',[250 1800],'qmax',.002,'qfloor',0,'side',1,'maxIter',10);
F=b1_forward_bi_p2c_fit(B,P,6,seed,[],struct('free',[]));
verifyLessThan(test,F.chi2,1e-10);
verifyLessThan(test,F.held_chi2n,1e-10);
verifyFalse(test,any(F.mask(:,1:2,:),'all'));
end
