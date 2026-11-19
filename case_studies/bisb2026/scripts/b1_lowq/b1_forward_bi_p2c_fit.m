function [F,D] = b1_forward_bi_p2c_fit(B,P,held,seed,L,control)
%B1_FORWARD_BI_P2C_FIT  Variable-projection fit for P2c on a matched operator.
% L is [] or a spectral column per unit material weight at Q_ref. Control
% optionally fixes material W and/or supplies a subset of nonlinear DL indices.
if nargin<5, L=[]; end
if nargin<6, control=struct(); end
t = tic;
ne = numel(B.cfg.nodes.E0);
ng = numel(B.cfg.nodes.G);
nj = B.intr.n_branch;
na = nj*numel(B.cfg.nodes.A);
if isfield(control,'design')
    D=control.design;
    limit=D.limit;
else
    [X,limit] = b1_forward_bi_p2c_gain(B,P.gain);
end
nk = numel(limit);
nDL = nj*(ne+ng);
u0 = [seed.u(1:nDL),zeros(1,nk)];
if isfield(seed,'gain') && numel(seed.gain)==nk, u0(nDL+1:end) = seed.gain; end
lb = [repelem(log(B.cfg.bounds_E0(1:nj,1).'.^2),ne),repmat(log(B.cfg.bounds_G(1)),1,nj*ng),-limit];
ub = [repelem(log(B.cfg.bounds_E0(1:nj,2).'.^2),ne),repmat(log(B.cfg.bounds_G(2)),1,nj*ng),limit];
u0 = min(max(u0,lb),ub);
mask = B.mask & B.E>=P.window(1) & B.E<=P.window(2) & abs(B.q)<=P.qmax+1e-10 & abs(B.q)>=P.qfloor-1e-10;
if isfield(P,'side'), mask = mask & sign(B.q)==P.side; end
train=mask; train(:,:,held)=false;
if ~isfield(control,'design')
    D = local_design(B,P,train,held);
    D.nDL = nDL;
    D.X = reshape(X,numel(B.Y),nk);
    D.limit = limit;
end
D.L = L;
D.fixedW = NaN;
if isfield(control,'W'), D.fixedW = control.W; end
D.fixedRatio = NaN;
if isfield(control,'ratio'), D.fixedRatio = control.ratio; end
free = 1:numel(u0);
if isfield(control,'free'), free = control.free; end
if isfield(P,'fixedGain') && ~isempty(P.fixedGain)
    assert(numel(P.fixedGain)==nk);
    u0(nDL+1:end) = P.fixedGain;
    free = free(free<=nDL);
end
D.noDLchange = ~any(free<=nDL);
if D.noDLchange
    if isfield(control,'cols0'), D.cols0=control.cols0; else, D.cols0=local_cols(u0,B,nDL); end
end
f = @(u) local_res(local_insert(u0,free,u),B,D);
nobs = nnz(train);
if isempty(free)
    u=u0; flag=1; it=0; J=zeros(nobs,0);
else
    maxIter=P.maxIter;
    if isfield(control,'maxIter'), maxIter=control.maxIter; end
    opt=optimoptions('lsqnonlin','Display','off','MaxIterations',maxIter,'MaxFunctionEvaluations',8000, ...
        'FunctionTolerance',1e-9,'StepTolerance',1e-7,'TypicalX',max(abs(u0(free)),.1));
    [v,~,~,flag,info,~,J]=lsqnonlin(f,u0(free),lb(free),ub(free),opt);
    u=local_insert(u0,free,v); it=info.iterations;
end
lite = isfield(control,'lite') && control.lite;
if lite
    [r,a]=local_res(u,B,D);
else
    [r,a,extra,g,m] = local_res(u,B,D);
end
F = struct('u',u(1:nDL),'gain',u(nDL+1:end),'a',a(1:na),'chi2',sum(r(1:nobs).^2), ...
    'objective',sum(r.^2),'held',held,'flag',flag,'iterations',it,'elapsed',toc(t));
th = b1_p2b_decode(F.u,B.cfg.nodes,nj);
F.E0 = sqrt(th.E0sq); F.G = th.G;
F.W = 0;
if ~isempty(L)
    F.W = D.fixedW;
    if ~isnan(D.fixedRatio)
        F.W = D.fixedRatio*(local_low_weight(th,B)*F.a(1:na/nj));
    elseif isnan(F.W)
        F.W = a(na+1);
    end
end
dof = nobs-numel(free)-numel(a)-D.nSigned;
if strcmp(D.stage,'low'), dof=dof-size(D.bg.A,2); end
F.chi2dof = F.chi2/dof;
F.dof = dof;
if strcmp(D.stage,'low')
    F.Wlow = local_low_weight(th,B)*F.a(1:na/nj);
    F.relativeW = F.W/max(F.Wlow,realmin);
end
if lite, return; end
V = zeros(numel(u),numel(u));
if ~isempty(free)
    V(free,free) = pinv(full(J.'*J))*max(1,F.chi2dof);
end
F.sE0 = F.E0.*reshape(sqrt(max(diag(V(1:nj*ne,1:nj*ne)),0)),ne,nj).'/2;
ig = nj*ne+(1:nj*ng);
F.sG = F.G.*reshape(sqrt(max(diag(V(ig,ig)),0)),ng,nj).';
F.sGain = sqrt(max(diag(V(nDL+1:end,nDL+1:end)),0)).';
F.n_train = nobs;
F.n_held = nnz(mask & ~train);
F.parts = zeros([size(B.Y) nj]);
for j=1:nj
    idx = (j-1)*na/nj+(1:na/nj);
    F.parts(:,:,:,j)=reshape(m(:,idx)*a(idx),size(B.Y));
end
F.mat = zeros(size(B.Y));
if ~isempty(L), F.mat = reshape(L(:).*g*F.W,size(B.Y)); end
F.tailCoef = extra.tailCoef; F.phCoef = extra.phCoef;
F.tail = extra.tail; F.phonon = extra.phonon;
F.mu = sum(F.parts,4)+F.mat+F.tail+F.phonon;
F.mask = mask;
Z = (F.mu-B.Y).^2./B.S2;
mv = train & B.E<=300;
F.train_valley = mean(Z(mv));
Z = Z(:,:,held);
mheld=mask(:,:,held); F.held_chi2n=mean(Z(mheld));
zwin=B.E<=300;
mvalley=mheld(zwin,:,:); zvalley=Z(zwin,:,:);
F.held_valley=mean(zvalley(mvalley));
end

function u=local_insert(u,idx,v)
u(idx)=v;
end

function A=local_cols(u,B,nDL)
th = b1_p2b_decode(u(1:nDL),B.cfg.nodes,B.intr.n_branch);
cols=qe_lowq_joint_forward('eval',B.cache,th,B.intr);
A=reshape(permute(cols,[1 3 4 2]),[],size(cols,2));
end

function w=local_low_weight(th,B)
E=4:4:2400;
q=B.intr.Q_ref;
C=qe_lowq_joint_forward('C',q,E,th,B.intr);
K=q./(q^2+(E*B.intr.qE_per_meV).^2).^2;
undo=-expm1(-E/B.intr.kT)./K;
n=numel(B.intr.Qn_A);
loss=C(:,:,1:n).*undo;
w=reshape(trapz([0 E],cat(2,zeros(1,1,n),loss),2),1,[])/1000;
end

function D=local_design(B,P,train,held)
[nE,nq,nb]=size(B.Y);
D = struct('stage',B.stage,'train',train,'y',B.Y(train),'w',1./sqrt(B.S2(train)), ...
    'held',held,'scale',B.exposure(:)/sum(B.exposure(setdiff(1:nb,held))), ...
    'lambda',B.cfg.lambda_curv,'nSigned',0);
if strcmp(B.stage,'odd')
    D.tail = (max(B.E,B.cfg.holdout(1))/B.cfg.tail.E_ref).^(-B.cfg.tail.p);
    D.tail(B.E<B.cfg.holdout(1))=0;
    idx = repmat(reshape(1:nq*nb,1,nq,nb),nE,1,1);
    D.s = idx(train);
    D.nSigned = numel(unique(D.s));
    T = repmat(D.tail,1,nq,nb);
    D.tw = T(train).*D.w;
    D.tt = accumarray(D.s,D.tw.^2,[nq*nb 1]);
    D.tt(D.tt==0)=Inf;
else
    D.bg = b1_forward_bi_p2c_background(B,P,held,train);
end
end

function [r,a,extra,g,m] = local_res(u,B,D)
if D.noDLchange, m=D.cols0; else, m=local_cols(u,B,D.nDL); end
g = exp(D.X*u(D.nDL+1:end).');
m = m.*g;
A = m(D.train(:),:);
y = D.y;
th = b1_p2b_decode(u(1:D.nDL),B.cfg.nodes,B.intr.n_branch);
if ~isempty(D.L)
    line = g(D.train(:)).*D.L(D.train);
    if ~isnan(D.fixedRatio)
        wb=local_low_weight(th,B);
        A(:,1:numel(wb))=A(:,1:numel(wb))+line*(D.fixedRatio*wb);
    elseif isnan(D.fixedW)
        A = [A line];
    else
        y = y-D.fixedW*line;
    end
end
Aw = A.*D.w; yw = y.*D.w;
pen = D.lambda*reshape(diff(log(th.E0sq),2,2),[],1);
if strcmp(D.stage,'odd')
    nt = numel(D.tt);
    ty = accumarray(D.s,D.tw.*yw,[nt 1]);
    tA = zeros(nt,size(A,2));
    for j=1:size(A,2), tA(:,j)=accumarray(D.s,D.tw.*Aw(:,j),[nt 1]); end
    Ap = Aw-D.tw.*tA(D.s,:)./D.tt(D.s);
    yp = yw-D.tw.*ty(D.s)./D.tt(D.s);
    a = local_nnls(Ap,yp,false);
    r = [Ap*a-yp;pen];
    if nargout>2
        coef = reshape((ty-tA*a)./D.tt, size(B.Y,2),size(B.Y,3));
        b = setdiff(1:size(B.Y,3),D.held);
        coef(:,D.held) = sum(coef(:,b),2)*D.scale(D.held);
        extra = struct('tailCoef',coef,'phCoef',[],'tail',D.tail.*reshape(coef,1,size(coef,1),[]), ...
            'phonon',zeros(size(B.Y)));
    end
else
    nLin = size(A,2);
    G = D.bg;
    AA = [Aw G.A;zeros(size(G.R,1),nLin) G.R];
    yy = [yw; G.target];
    allA = local_nnls(AA,yy,true);
    a = allA(1:nLin);
    r = [AA*allA-yy;pen];
    if nargout>2
        extra = b1_forward_bi_p2c_background(B,G,allA(nLin+1:end),D.held);
    end
end
end

function a=local_nnls(A,y,compress)
sc = max(abs(A),[],1);
sc(sc==0)=1;
A=A./sc;
if compress
    [R,p]=chol(A.'*A);
    if p==0 && rcond(R)>1e-5
        y=R.'\(A.'*y);
        A=R;
    end
end
a = lsqnonneg(A,y)./sc(:);
end
