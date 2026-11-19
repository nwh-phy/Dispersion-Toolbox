function tests=test_b1_forward_bi_p2c_search
tests=functiontests(localfunctions);
end

function testMaterialRecoveryAndStrictHoldout(test)
[B,P,seed,L]=local_fixture;
[M,T,U]=b1_forward_bi_p2c_search(B,P,6,seed,true);
verifyEqual(test,M{2}.shape,[200 25]);
verifyEqual(test,M{2}.W,3e-3,'RelTol',.01);
verifyGreaterThan(test,M{1}.chi2,100*M{2}.chi2);
verifyTrue(test,any(T.refit));
verifyGreaterThan(test,max(U.upper(strcmp(U.metric,'W') & strcmp(U.family,'interband'))),M{2}.W);
B.Y(B.E>=P.window(1),:,6)=B.Y(B.E>=P.window(1),:,6)+15;
S=b1_forward_bi_p2c_fit(B,P,6,M{2},L,struct('free',[]));
verifyEqual(test,S.mu,M{2}.mu,'RelTol',1e-9,'AbsTol',1e-7);
verifyGreaterThan(test,S.held_chi2n,M{2}.held_chi2n);
end

function testProfileBracketsInSmallWeightUnits(test)
[B,P,seed]=local_fixture;
B.cache.H=B.cache.H*1e12;
P.nesting_check=true;
[M,~,U]=b1_forward_bi_p2c_search(B,P,6,seed,true);
verifyEqual(test,M{2}.W,3e-15,'RelTol',.01);
q=strcmp(U.metric,'W') & strcmp(U.family,'interband') & U.refit;
verifyTrue(test,all(U.closed(q)));
verifyLessThan(test,max(U.upper(q)),1e-13);
verifyLessThan(test,max(U.bracket_delta(q)),.15);
end

function [B,P,seed,L] = local_fixture
q=[-.003 -.002 -.001 .001 .002 .003];
Q=[0;.001;.002;.003];
W=zeros(6,1,4);
for i=1:numel(q), W(i,1,round(1+abs(q(i))*1000))=1e-6; end
nb=6;
ops=repmat({struct('W',W,'Qm',Q)},1,nb);
C=qe_lowq_joint_forward('build',ops,repmat({[0 1 0]},1,nb), ...
    struct('omega',-400:4:1800,'eps',-4:4:4,'variant','slice'));
E=C.E_out(:);
nodes=struct('E0',[.001 .002 .003],'G',[.001 .003],'A',[.001 .003]);
I=qe_lowq_joint_forward('intrinsic','Qn_E0',nodes.E0,'Qn_G',nodes.G,'Qn_A',nodes.A,'n_branch',1);
th=struct('E0sq',[450 560 640].^2,'G',[210 350]);
a=qe_lowq_joint_forward('eval',C,th,I);
y=reshape(sum(a.*reshape([.4 .6]*1e8,1,[],1,1),2),numel(E),6,6);
cfg=struct('nodes',nodes,'bounds_E0',[150 950],'bounds_G',[20 2400],'lambda_curv',0);
k=struct('c',0,'A',[1 0],'sL',[11 11],'sR',[11 11],'mL',[2 2],'mR',[2 2]);
R=struct('centre',zeros(nb,6),'area',1000*ones(nb,6),'th',{repmat({k},nb,6)});
z=repmat(30*(100./max(abs(E),30)).^3,1,6,nb);
B=struct('E',E,'Y',y,'S2',max(y,1),'q',q,'cfg',cfg,'intr',I, ...
    'cache',C,'mask',true(size(y)),'exposure',ones(nb,1),'stage','low','response',R, ...
    'cube',struct('E',E,'q',q,'b6',struct('J',z)));
P=struct('id','synthetic','gain','none','window',[200 1800],'qmax',.003,'qfloor',0, ...
    'maxIter',50,'fullB1',false,'briefGrid',true,'tail','mirror','phonon',false,'gainBands',[70 180;85 240;100 320]);
m=repmat(E>=200 & E<=1800,1,6,nb); m(:,:,6)=false;
D=b1_forward_bi_p2c_background(B,P,6,m);
L=b1_forward_bi_p2c_material(B,'interband',[200 25],200);
B.Y=y+D.tau+3e-3*L;
B.S2=max(B.Y,1);
seed=struct('u',[log(th.E0sq) log(th.G)]);
end
