function fits = qe_compare_component_models(E, Y, options)
% Same-window explicit n=1/2 bounded multistart pilot; no deletion or repair.
% Uses peak_models definitions. Unweighted LS is a descriptive objective,
% not a calibrated likelihood for corrected/interpolated detector data.
arguments
 E (:,1) double
 Y (:,1) double
 options.energy_window (1,2) double = [300 1800]
 options.peak_model char = 'lorentz_symmetric'
 options.n_starts (1,1) double {mustBePositive,mustBeInteger} = 24
 options.seed (1,1) double = 20260912
 options.baseline_mode char {mustBeMember(options.baseline_mode,{'power_law','power_law_plus_nonnegative_constant'})} = 'power_law'
 options.start_policy char {mustBeMember(options.start_policy,{'legacy','independent'})} = 'legacy'
 options.extra_starts (1,1) double {mustBeNonnegative,mustBeInteger} = 0
 options.h0_n_starts (1,1) double = NaN
 options.max_iterations (1,1) double {mustBeNonnegative,mustBeInteger} = 600
end
assert(ismember(options.peak_model,{'lorentz','lorentz_symmetric'}),'Unsupported pilot model');
assert(numel(E)==numel(Y)&&all(diff(E)>0),'Invalid energy axis');
mask=E>=options.energy_window(1)&E<=options.energy_window(2)&isfinite(Y);
E=E(mask); Y=Y(mask); assert(numel(E)>=12,'Insufficient common-window samples');
model=peak_models(options.peak_model); scale=max(abs(Y)); if scale==0, scale=1; end
ampunit=1000; if strcmp(options.peak_model,'lorentz'), ampunit=1e6; end
stream=RandStream('mt19937ar','Seed',options.seed);
opt=optimoptions('lsqnonlin','Display','off','MaxIterations',options.max_iterations, ...
 'MaxFunctionEvaluations',12000,'FunctionTolerance',1e-10,'StepTolerance',1e-10);
fits=struct([]);
for n=1:2
 if strcmp(options.start_policy,'independent'), stream=RandStream('mt19937ar','Seed',options.seed+1009*n); end
 lb=[0 0 repmat([min(E)/1000 max(median(diff(E)),1)/1000 0],1,n)];
 ub=[Inf 6 repmat([max(E)/1000 5 Inf],1,n)];
 if strcmp(options.baseline_mode,'power_law_plus_nonnegative_constant'), lb(end+1)=0; ub(end+1)=Inf; end
 candidates=struct([]);
 count=options.n_starts;
 if n==1 && isfinite(options.h0_n_starts), count=options.h0_n_starts; end
 if n==2, count=count+options.extra_starts; end
 for k=1:count
  c=sort(min(E)/1000+(max(E)-min(E))/1000*(.05+.9*rand(stream,1,n)));
  if k<=4
   c=linspace(min(E)/1000+.2,max(E)/1000-.2,n+2); c=c(2:end-1);
  end
  widths=[.08 .25 .65 1.5]; w=widths(1+mod(k-1,4))*(.6+.8*rand(stream,1,n));
  p0=[max(.001,min(Y/scale))*(.3+.7*rand(stream)) .2+3*rand(stream)];
  for j=1:n, p0=[p0 c(j) w(j) .1+rand(stream)*2]; end %#ok<AGROW>
  if n==2 && k>options.n_starts
   p0(4)=.04+.08*rand(stream); p0(7)=.6+2*rand(stream);
   if mod(k,2)==0, p0([4 7])=p0([7 4]); end
   p0(8)=1e-5+.03*rand(stream);
  end
  if numel(lb)>2+3*n, p0(end+1)=.01; end
  p0=max(lb+1e-8,min(ub-1e-8,p0));
  row=struct('start',k,'p0',p0,'p',nan(size(p0)),'exitflag',NaN, ...
   'objective',Inf,'boundary',false(size(p0)),'jacobian_condition',NaN,'message','', ...
   'firstorderopt',NaN,'iterations',NaN,'funcCount',NaN,'lambda',struct(), ...
   'jacobian',[],'singular_values',[],'column_scaled_singular_values',[], ...
   'candidate_type',options.start_policy);
  if k>options.n_starts, row.candidate_type='asymmetric_weak_extra'; end
  try
   assert(all(isfinite(p0))&&all(p0>=lb)&all(p0<=ub),'Invalid initial parameters');
   [p,obj,~,flag,output,lambda,J]=lsqnonlin(@(p) prediction(p)-Y/scale,p0,lb,ub,opt);
   row.p=p; row.exitflag=flag; row.objective=obj; row.message=output.message;
   row.boundary=abs(p-lb)<1e-5.*max(1,abs(lb)) | abs(p-ub)<1e-5.*max(1,abs(ub));
   row.boundary(~isfinite(ub))=abs(p(~isfinite(ub))-lb(~isfinite(ub)))<1e-5;
   row.jacobian_condition=cond(full(J));
   row.firstorderopt=output.firstorderopt; row.iterations=output.iterations; row.funcCount=output.funcCount;
   row.lambda=lambda; row.jacobian=full(J); row.singular_values=svd(full(J));
   row.column_scaled_singular_values=svd(full(J)./max(vecnorm(full(J)),eps));
  catch ME
   row.message=[ME.identifier ': ' ME.message];
  end
  if k==1, candidates=row; else, candidates(k)=row; end
 end
 good=find([candidates.exitflag]>0 & isfinite([candidates.objective]));
 f=struct('n_components',n,'peak_model',options.peak_model,'energy_meV',E, ...
  'observed',Y,'scale',scale,'ampunit',ampunit,'lb',lb,'ub',ub,'candidates',candidates, ...
  'success',false,'selected_start',NaN,'parameters',nan(n,3),'background',nan(size(E)), ...
  'components',nan(numel(E),n),'prediction',nan(size(E)),'residual',nan(size(E)), ...
  'sse',NaN,'normalized_sse',NaN,'numerical_status','failed', ...
  'scientific_status','unresolved_or_invalid','uncertainty_status','not_run_P4_P5', ...
  'fwhm_diagnostic',struct([]),'baseline_mode',options.baseline_mode,'effective_options',options, ...
  'component_order',[],'raw_to_order',[],'witness',struct('available',false), ...
  'identifiability_status','not_assessed','origin_status','not_assessed');
 if ~isempty(good)
  [~,ki]=min([candidates(good).objective]); best=good(ki); p=candidates(best).p;
  [~,bg,components]=prediction(p);
  pars=reshape(p(3:2+3*n),3,n).'; [~,order]=sort(pars(:,1)); pars=pars(order,:);
  f.component_order=order; f.raw_to_order=zeros(1,n); f.raw_to_order(order)=1:n;
  f.parameters=pars.*[1000 1000 ampunit*scale];
  f.background=bg*scale; f.components=components(:,order)*scale;
  f.prediction=f.background+sum(f.components,2); f.residual=Y-f.prediction;
  f.sse=sum(f.residual.^2); f.normalized_sse=f.sse/scale^2;
  assert(abs(f.normalized_sse-candidates(best).objective)<1e-7*max(1,f.normalized_sse),'Objective contract');
  f.success=true; f.selected_start=best; f.numerical_status='converged';
  if any(candidates(best).boundary), f.numerical_status='boundary'; end
  if n==2 && diff(f.parameters(:,1))<median(diff(E)), f.numerical_status='collapse_sampling_scale'; end
  for j=1:n
   measurement=measure_peak_fwhm(E,f.components(:,j));
   if j==1, f.fwhm_diagnostic=measurement; else, f.fwhm_diagnostic(j)=measurement; end
  end
 end
 if n==2 && fits(1).success
  h0=fits(1); p1=h0.candidates(h0.selected_start).p;
  pw=[p1(1:5) mean(E)/1000 .25 0];
  if numel(p1)>5, pw(end+1)=p1(end); end
  yw=prediction(pw)*scale; qw=sum((yw-Y).^2)/scale^2;
  f.witness=struct('available',true,'p',pw,'prediction',yw,'objective',qw, ...
   'h0_objective',h0.normalized_sse,'candidate_type','feasible_zero_amplitude_not_optimized', ...
   'best_feasible_objective',min(qw,f.normalized_sse), ...
   'optimized_nesting_violation',~f.success || f.normalized_sse>h0.normalized_sse+1e-7*max(1,h0.normalized_sse));
 end
 if n==1, fits=f; else, fits(n)=f; end
end
 function [total,bg,comp]=prediction(p)
  bg=p(1)*(E/1000).^(-p(2)); comp=zeros(numel(E),n);
  if numel(p)>2+3*n, bg=bg+p(end); end
  for jj=1:size(comp,2)
   ii=3+(jj-1)*3;
   comp(:,jj)=model.model_fn(p(ii)*1000,p(ii+1)*1000,p(ii+2)*ampunit,E);
  end
  total=bg+sum(comp,2);
 end
end
