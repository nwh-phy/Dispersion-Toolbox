function G = b1_forward_bi_p2c_background(B,P,held,train)
%B1_FORWARD_BI_P2C_BACKGROUND  Detector tail and 130 meV block templates.
% Evaluation form: background(B,design,coefficients,heldBlock).
if nargin==4 && isfield(P,'tau')
    G = local_prediction(B,P,held,train);
    return
end
E=B.E;
[nE,ni,nb]=size(B.Y);
dE=median(diff(B.cube.E));
[tf,idx]=ismembertol(B.q,B.cube.q,1e-10,'DataScale',1);
assert(all(tf));
T=zeros(nE,ni,nb); Ph=T;
J=B.cube.b6.J;
gains=P.gainBands;
pow=find(strcmp(P.tail,{'power1','power2','power3'}));
zero=P.phonon==false;
prior=zeros(1,nb);
priorDev=ones(1,nb);
for b=1:nb
    pulse=zeros(1,ni); sigma=pulse; alt=pulse;
    for i=1:ni
        j=idx(i);
        centre=B.response.centre(b,j);
        x=B.cube.E(:)-centre;
        y=double(J(:,j,b));
        [T(:,i,b),mirror]=local_tail(E-centre,x,y,P.tail,pow,gains);
        if zero, continue; end
        th=B.response.th{b,j};
        a=B.response.area(b,j);
        % Coefficient is pulse area in units of 10^-5 times local elastic area.
        Ph(:,i,b)=qe_elastic_line_model('kernel',E-130,th,60)*dE*a*1e-5;
        if b==held, continue; end     % no held loss-side phonon calibration
        [pulse(i),sigma(i),alt(i)]=local_phonon(x,y,mirror,th,dE,a*1e-5);
    end
    if zero || b==held, continue; end
    in=any(train(:,:,b),1) & isfinite(sigma) & sigma>0;
    if ~any(in), continue; end
    % Add a model-spread floor when pooling channel estimates. The twice
    % sideband/model spread is a regularizing scale, not a confidence interval.
    s=hypot(sigma(in),abs(alt(in)));
    w=1./max(s,0.3).^2;
    pm=sum(pulse(in).*w)/sum(w);
    prior(b)=max(0,pm);
    dev=sqrt(sum(w.*(pulse(in)-pm).^2)/sum(w));
    priorDev(b)=max([.3,2/sqrt(sum(w)),2*dev,abs(pm)*.5]);
end
blocks=setdiff(1:nb,held);
ntr=numel(blocks);
nt=ntr*(1+~zero);
R=zeros(nt);
target=zeros(nt,1);
A=zeros(nnz(train),nt);
maskB=reshape(1:nb,1,1,nb)==reshape(blocks,1,1,1,ntr);
S=B.S2(train);
for k=1:ntr
    v=T.*maskB(:,:,:,k);
    A(:,k)=v(train)./sqrt(S);
    R(k,k)=1/.75;
    target(k)=R(k,k);
    if ~zero
        v=Ph.*maskB(:,:,:,k);
        A(:,k+ntr)=v(train)./sqrt(S);
        R(k+ntr,k+ntr)=1/priorDev(blocks(k));
        target(k+ntr)=prior(blocks(k))*R(k+ntr,k+ntr);
    end
end
G=struct('A',A,'R',R,'target',target,'tau',T,'ph',Ph,'zero',zero, ...
    'blocks',blocks,'ph_prior',prior,'ph_scale',priorDev);
end

function [t,yGain] = local_tail(E,x,y,name,pow,gains)
% Reflection uses the measured gain side around each slice centre. At high
% energy match its integrated 280--360 meV level onto an E^-4 continuation.
% Smooth only the reflected noise. The 68--84 meV anchor is disjoint from
% all fitted windows and the phonon band.
good=isfinite(y);
xg=x(good); jg=y(good);
yGain=interp1(xg,jg,-x,'linear',0);
z=max(E,4);
if isempty(pow)
    ys=movmean(jg,5,1);
    m=(-xg>=280 & -xg<=360);
    ref=median(max(ys(m),0).*(-xg(m)/320).^4);
    if isempty(ref) || ~isfinite(ref), ref=0; end
    mirror=interp1(xg,ys,-z,'linear',0);
    continuation=ref*(z/320).^-4;
    blend=min(max((z-280)/80,0),1);
    t=(1-blend).*max(mirror,0)+blend.*continuation;
    if strcmp(name,'anchor')
        a=x>=68 & x<=84 & isfinite(y);
        t=t*max(sum(y(a))/max(sum(yGain(a)),eps),0);
    end
else
    [~,p,A]=qe_powerlaw_tail(-x,y,gains(pow,:));
    assert(isfinite(p) && isfinite(A) && A>=0,'Invalid gain-side power law.');
    blend=min(max((z-320)/80,0),1);
    t=(1-blend).*A.*z.^(-p)+blend.*A*320^(-p).*(z/320).^-4;
end
t=t.*min(max((z-60)/12,0),1);
end

function [p,s,sdiff]=local_phonon(x,y,yg,th,dE,scale)
if scale<=0, p=0; s=Inf; sdiff=0; return; end
u=x>=112 & x<=148 & isfinite(y);
side=[80 108;156 188];
m=(x>=side(1,1) & x<=side(1,2)) | (x>=side(2,1) & x<=side(2,2));
alt=(x>=88 & x<=112) | (x>=164 & x<=204);
p0=polyfit(log(x(m)),log(max(y(m),1)),1);
p1=polyfit(log(x(alt)),log(max(y(alt),1)),1);
p2=polyfit(log(x(m)),log(max(yg(m),1)),1);
th.c=0;
g=qe_elastic_line_model('kernel',x(u)-130,th,60)*dE*scale;
w=1./max(y(u),1);
d=sum(w.*g.^2);
f=@(z,t)sum(w.*g.*(z(u)-exp(polyval(t,log(x(u))))))/d;
p=f(y,p0);
s=sqrt(1/d);
sdiff=max(abs(p-f(y,p1)),abs(f(yg,p2)));
end

function G=local_prediction(B,D,c,held)
nb=size(B.Y,3);
n=numel(D.blocks);
tail=zeros(1,nb); ph=tail;
tail(D.blocks)=c(1:n);
% Each tau already carries the held gain-side amplitude for its own block.
% Tail coefficients are correction factors. Phonon coefficients are
% fractions of contemporaneous elastic area. Transfer an exposure-weighted rate.
w=B.exposure(D.blocks);
tail(held)=sum(tail(D.blocks).*w(:).')/sum(w);
if ~D.zero
    ph(D.blocks)=c(n+1:end);
    ph(held)=sum(ph(D.blocks).*w(:).')/sum(w);
end
G=struct('tailCoef',tail,'phCoef',ph,'tail',D.tau.*reshape(tail,1,1,[]), ...
    'phonon',D.ph.*reshape(ph,1,1,[]));
end
