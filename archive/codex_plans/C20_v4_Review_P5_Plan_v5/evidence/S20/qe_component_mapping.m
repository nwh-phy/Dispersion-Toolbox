function m = qe_component_mapping(p,lb,ub,n,scale,ampunit,peak_model)
% Optimizer slots -> energy-ordered labels, with original v2 tolerance.
p=p(:).'; lb=lb(:).'; ub=ub(:).';
assert(numel(p)==numel(lb)&&numel(p)==numel(ub));
pars=reshape(p(3:2+3*n),3,n).';
[~,order]=sort(pars(:,1)); inverse=zeros(1,n); inverse(order)=1:n;
native_scale=[scale 1 repmat([1000 1000 ampunit*scale],1,n)];
names=["background_B","background_r",repmat(["E0","width","A"],1,n)];
units=["counts","dimensionless",repmat(["meV","meV","counts*meV"],1,n)];
if strcmp(peak_model,'lorentz'), units(5:3:2+3*n)="counts*meV^2"; end
ordered=zeros(size(p)); raw_component=zeros(size(p));
for j=1:n
 ii=2+3*(j-1)+(1:3); ordered(ii)=inverse(j); raw_component(ii)=j;
end
if numel(p)>2+3*n
 native_scale(end+1)=scale; names(end+1)="background_C"; units(end+1)="counts";
end
tol_lo=1e-5*max(1,abs(lb)); tol_hi=1e-5*max(1,abs(ub));
tol_lo(~isfinite(ub))=1e-5;
lower=isfinite(p)&abs(p-lb)<tol_lo;
upper=isfinite(p)&isfinite(ub)&abs(p-ub)<tol_hi;
types=repmat("interior",size(p));
types(lower)=names(lower)+"_lower"; types(upper)=names(upper)+"_upper";
types(~isfinite(p))="not_assessable";
m=struct('order',order,'raw_to_order',inverse,'raw_component',raw_component, ...
 'ordered_component',ordered,'parameters',pars(order,:).*[1000 1000 ampunit*scale], ...
 'native_scale',native_scale,'names',names,'units',units,'lower_hit',lower,'upper_hit',upper, ...
 'tol_lower',tol_lo,'tol_upper',tol_hi,'boundary_type',types,'finite',all(isfinite(p)), ...
 'equal_center',n==2 && abs(pars(1,1)-pars(end,1))<1e-8);
end
