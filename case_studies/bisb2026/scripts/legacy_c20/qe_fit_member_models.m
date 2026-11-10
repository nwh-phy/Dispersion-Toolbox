function fits = qe_fit_member_models(E,q,X,options)
arguments
 E (:,1) double
 q (1,:) double
 X double
 options.energy_window (1,2) double = [300 1800]
 options.peak_model char {mustBeMember(options.peak_model,{'lorentz','lorentz_symmetric'})} = 'lorentz_symmetric'
 options.n_starts (1,1) double {mustBePositive,mustBeInteger} = 12
 options.seed (1,1) double = 50912
 options.max_iterations (1,1) double {mustBeNonnegative,mustBeInteger} = 600
 options.initial_fits struct = struct([])
end
assert(size(X,1)==numel(E)&&size(X,2)==numel(q)&&numel(q)>=3&&all(diff(q)>0));
mask=E>=options.energy_window(1)&E<=options.energy_window(2)&all(isfinite(X),2);
E=E(mask); X=X(mask,:); assert(numel(E)>=12);
scale=max(abs(X),[],'all'); if scale==0, scale=1; end
Q=numel(q); ampunit=1000; if strcmp(options.peak_model,'lorentz'), ampunit=1e6; end
opts=optimoptions('lsqnonlin','Display','off','SpecifyObjectiveGradient',true,'MaxIterations',options.max_iterations, ...
 'MaxFunctionEvaluations',12000,'FunctionTolerance',1e-10,'StepTolerance',1e-10,'OptimalityTolerance',1e-7);
fits=struct([]);
for n=1:2
 stream=RandStream('mt19937ar','Seed',options.seed+1009*n);
 lb=[0 zeros(1,Q) repmat([min(E)/1000 min(E)/1000 max(median(diff(E)),1)/1000 zeros(1,Q)],1,n)];
 ub=[6 inf(1,Q) repmat([max(E)/1000 max(E)/1000 5 inf(1,Q)],1,n)]; candidates=struct([]);
 for start=1:options.n_starts
  centers=sort(.4+1.2*rand(stream,1,n)); widths=[.2 .8 1.5 3];
  p0=[rand(stream)*2 .03*ones(1,Q)];
  for j=1:n, p0=[p0 centers(j)+.15*(rand(stream,1,2)-.5) widths(1+mod(start-1,4)) .5*ones(1,Q)]; end %#ok<AGROW>
  if start==1 && ~isempty(options.initial_fits)&&options.initial_fits(n).success
   f=options.initial_fits(n); z=f.candidates(f.selected_start).p; p0(1)=z(2); p0(2:Q+1)=z(1)*f.scale/scale;
   for j=1:n
    base=1+Q+(j-1)*(3+Q); pp=f.parameters(j,:);
    p0(base+(1:3))=[pp(1)/1000 pp(1)/1000 pp(2)/1000]; p0(base+3+(1:Q))=pp(3)/ampunit/scale;
   end
  end
  p0=max(lb+1e-8,min(ub-1e-8,p0));
  c=struct('p0',p0,'p',nan(size(p0)),'objective',Inf,'exitflag',NaN,'output',struct(), ...
   'lambda',struct(),'J_singular_values',[],'scaled_J_singular_values',[],'boundary',false(size(p0)),'order',[],'message','');
  try
   [p,v,~,flag,output,lambda,J]=lsqnonlin(@objective,p0,lb,ub,opts);
   c.p=p; c.objective=v; c.exitflag=flag; c.output=output; c.lambda=lambda;
   c.J_singular_values=svd(full(J)); c.scaled_J_singular_values=svd(full(J)./max(vecnorm(full(J)),eps));
   c.boundary=abs(p-lb)<1e-5; c.boundary(isfinite(ub))=c.boundary(isfinite(ub))|abs(p(isfinite(ub))-ub(isfinite(ub)))<1e-5;
   center=zeros(1,n); for j=1:n, b=1+Q+(j-1)*(3+Q); center(j)=mean(p(b+(1:2))); end
   [~,c.order]=sort(center); c.message=output.message;
  catch ME, c.message=[ME.identifier ': ' ME.message]; end
  if isempty(candidates), candidates=c; else, candidates(end+1)=c; end %#ok<AGROW>
 end
 f=struct('E',E,'q',q,'observed',X,'scale',scale,'ampunit',ampunit,'n_components',n,'peak_model',options.peak_model, ...
  'effective_options',options,'lb',lb,'ub',ub,'candidates',candidates,'success',false,'selected_start',NaN,'p',nan(size(lb)), ...
  'prediction',nan(size(X)),'background',nan(size(X)),'components',nan(size(X,1),Q,n),'residual',nan(size(X)), ...
  'objective',NaN,'order',[],'native_centers',nan(Q,n),'native_widths',nan(1,n),'native_A',nan(Q,n), ...
  'slope_meV_A',nan(1,n),'witness',struct('available',false),'observation_count',numel(X), ...
  'objective_input','five member spectra once; derived bins excluded','identifiability','not_assessed');
 good=find([candidates.exitflag]>0&isfinite([candidates.objective]));
 if ~isempty(good)
  [~,k]=min([candidates(good).objective]); k=good(k); p=candidates(k).p;
  [Y,B,C]=qe_member_prediction(E,q,p,n,options.peak_model,ampunit);
  f.success=true; f.selected_start=k; f.p=p; f.order=candidates(k).order;
  f.prediction=Y*scale; f.background=B*scale; f.components=C(:,:,f.order)*scale; f.residual=X-f.prediction;
  f.objective=sum(f.residual.^2,'all')/scale^2;
  beta=(q-min(q))/(max(q)-min(q));
  for j=1:n
   raw=f.order(j); b=1+Q+(raw-1)*(3+Q);
   f.native_centers(:,j)=1000*((1-beta)*p(b+1)+beta*p(b+2));
   f.native_widths(j)=1000*p(b+3); f.native_A(:,j)=p(b+3+(1:Q))*ampunit*scale;
   f.slope_meV_A(j)=1000*(p(b+2)-p(b+1))/(max(q)-min(q));
  end
 end
 if n==2&&fits(1).success
  p=[fits(1).p mean(E)/1000 mean(E)/1000 .5 zeros(1,Q)];
  yw=qe_member_prediction(E,q,p,2,options.peak_model,ampunit)*scale; qw=sum((yw-X).^2,'all')/scale^2;
  assert(abs(qw-fits(1).objective)<1e-7*max(1,qw));
  f.witness=struct('available',true,'p',p,'prediction',yw,'objective',qw,'type','M1_embedded_all_second_amplitudes_zero', ...
   'optimized_violation',~f.success||f.objective>qw+1e-7*max(1,qw));
 end
 if isempty(fits), fits=f; else, fits(n)=f; end
end
 function [r,J]=objective(p)
  [Y,~,~,J]=qe_member_prediction(E,q,p,n,options.peak_model,ampunit); r=Y-X/scale; r=r(:);
 end
end
