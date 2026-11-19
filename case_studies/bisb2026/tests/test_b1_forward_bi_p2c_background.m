function tests=test_b1_forward_bi_p2c_background
tests=functiontests(localfunctions);
end

function testHeldGainAndIndependentLoss(test)
E=(-400:4:600).'; q=[-.002 -.001 .001 .002]; nb=6;
th=struct('c',0,'A',[1 0],'sL',[11 11],'sR',[11 11],'mL',[2 2],'mR',[2 2]);
R=struct('centre',zeros(nb,4),'area',10000*ones(nb,4),'th',{repmat({th},nb,4)});
J=repmat(120*(100./max(abs(E),30)).^3,1,4,nb);
B=struct('E',E,'Y',J,'S2',max(J,1),'q',q,'cube',struct('E',E,'q',q,'b6',struct('J',J)), ...
    'response',R,'exposure',ones(nb,1));
P=struct('tail','mirror','phonon',true,'gainBands',[70 180;85 240;100 320]);
m=repmat(E>=150,1,4,nb); m(:,:,6)=false;
D=b1_forward_bi_p2c_background(B,P,6,m);
verifyEqual(test,D.blocks,1:5);
verifyEqual(test,D.ph_prior(6),0);
B.cube.b6.J(E>=100,:,6)=B.cube.b6.J(E>=100,:,6)*10;
Z=b1_forward_bi_p2c_background(B,P,6,m);
verifyEqual(test,D,Z);
C=b1_forward_bi_p2c_background(B,D,[1:5 3:7].',6);
verifyEqual(test,C.tailCoef(6),3);
verifyEqual(test,C.phCoef(6),5);
verifyGreaterThan(test,max(C.phonon,[],'all'),0);
end
