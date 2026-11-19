function tests=test_b1_forward_bi_p2c_material
tests=functiontests(localfunctions);
end

function testUnitWeightAndDirectConvolution(test)
Q=[0;.001;.002];
op=struct('Qm',Q,'W',reshape([0 2 0],1,1,3));
c=qe_lowq_joint_forward('build',{op},{[.2 .6 .2]}, ...
    struct('omega',-400:4:800,'eps',-4:4:4,'variant','slice'));
B=struct('cache',c,'E',c.E_out(:),'intr',qe_lowq_joint_forward('intrinsic'));
G=80;
[v,f,n]=b1_forward_bi_p2c_material(B,'drude',G,150);
fun=@(x) G*x./(x.^2+G^2)./(1+(x/1100).^2).^2;
expect=integral(fun,150,300)*1e-3;
verifyEqual(test,n,expect,'RelTol',1e-12);
verifyEqual(test,f,integral(fun,250,300)*1e-3/expect,'RelTol',1e-12);
ww=c.omega; a=max(abs(ww),1);
b=ww./(-expm1(-ww/B.intr.kT)); b(ww==0)=B.intr.kT;
kernel=2*.001./(.001^2+(a*B.intr.qE_per_meV).^2).^2;
z=conv(kernel.*b.*fun(a)./(n*a),[.2 .6 .2]);
verifyEqual(test,v,z(:),'RelTol',1e-7,'AbsTol',1e-4);
[~,~,n270]=b1_forward_bi_p2c_material(B,'drude',G,270);
verifyEqual(test,n270,integral(fun,270,300)*1e-3,'RelTol',1e-12);
end
