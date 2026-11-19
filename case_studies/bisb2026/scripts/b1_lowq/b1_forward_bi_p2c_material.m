function [L,frac,Wshape] = b1_forward_bi_p2c_material(B,family,pars,Elo)
%B1_FORWARD_BI_P2C_MATERIAL  Forward a unit valley-weight loss spectrum.
% The positive L at Q_ref integrates to one on Elo--300 meV when Elo<300.
% With Elo=300 the normalization remains 250--300 for the baseline control.
% frac is the corresponding weight on the common 250--300 meV interval.
a=Elo;
if Elo>=300, a=250; end
f=@(e) local_shape(e,family,pars);
Wshape=integral(f,a,300)/1000;
assert(isfinite(Wshape) && Wshape>0);
frac=integral(f,250,300)/1000/Wshape;
w=B.cache.omega;
aw=max(abs(w),abs(w(2)-w(1))/4);
g=w./(-expm1(-w/B.intr.kT)); g(w==0)=B.intr.kT;
q=B.cache.Qm(:);
qE=aw*B.intr.qE_per_meV;
K=q./(q.^2+qE.^2).^2;
C=K.*(q/B.intr.Q_ref).*((f(aw)./aw).*g/Wshape);
Lf=fft(C,B.cache.N,2);
z=pagemtimes(B.cache.H,reshape(Lf(:,1:B.cache.Nf),numel(q),1,[]));
z=cat(3,z,conj(z(:,:,B.cache.Nf-1:-1:2)));
y=real(ifft(z,B.cache.N,3));
L=permute(reshape(y(:,:,1:numel(B.E)),B.cache.n_exit,B.cache.n_block,numel(B.E)),[3 1 2]);
end

function y=local_shape(e,family,p)
% Both tails have a finite optical first moment and the same leading
% 2D q scaling (v_2D chi'' with chi'' proportional to q^2).
switch family
    case 'interband'
        y=(-expm1(-e/60))./(1+exp((p(1)-e)/p(2)))./(1+(e/1100).^4);
    case 'drude'
        y=(p(1)*e./(e.^2+p(1)^2))./(1+(e/1100).^2).^2;
    otherwise
        error('Unknown low-energy material family %s.',family);
end
end
