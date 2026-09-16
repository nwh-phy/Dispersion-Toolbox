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
end
assert(ismember(options.peak_model,{'lorentz','lorentz_symmetric'}),'Unsupported pilot model');
assert(numel(E)==numel(Y)&&all(diff(E)>0),'Invalid energy axis');
mask=E>=options.energy_window(1)&E<=options.energy_window(2)&isfinite(Y);
E=E(mask); Y=Y(mask); assert(numel(E)>=12,'Insufficient common-window samples');
model=peak_models(options.peak_model); scale=max(abs(Y)); if scale==0, scale=1; end
ampunit=1000; if strcmp(options.peak_model,'lorentz'), ampunit=1e6; end
stream=RandStream('mt19937ar','Seed',options.seed);
opt=optimoptions('lsqnonlin','Display','off','MaxIterations',600, ...
 'MaxFunctionEvaluations',12000,'FunctionTolerance',1e-10,'StepTolerance',1e-10);
fits=struct([]);
for n=1:2
 lb=[0 0 repmat([min(E)/1000 max(median(diff(E)),1)/1000 0],1,n)];
 ub=[Inf 6 repmat([max(E)/1000 5 Inf],1,n)];
 candidates=struct([]);
 for k=1:options.n_starts
  c=sort(min(E)/1000+(max(E)-min(E))/1000*(.05+.9*rand(stream,1,n)));
  if k<=4
   c=linspace(min(E)/1000+.2,max(E)/1000-.2,n+2); c=c(2:end-1);
  end
  widths=[.08 .25 .65 1.5]; w=widths(1+mod(k-1,4))*(.6+.8*rand(stream,1,n));
  p0=[max(.001,min(Y/scale))*(.3+.7*rand(stream)) .2+3*rand(stream)];
  for j=1:n, p0=[p0 c(j) w(j) .1+rand(stream)*2]; end %#ok<AGROW>
  p0=max(lb+1e-8,min(ub-1e-8,p0));
  row=struct('start',k,'p0',p0,'p',nan(size(p0)),'exitflag',NaN, ...
   'objective',Inf,'boundary',false(size(p0)),'jacobian_condition',NaN,'message','');
  try
   assert(all(isfinite(p0))&&all(p0>=lb)&all(p0<=ub),'Invalid initial parameters');
   [p,obj,~,flag,output,~,J]=lsqnonlin(@(p) prediction(p)-Y/scale,p0,lb,ub,opt);
   row.p=p; row.exitflag=flag; row.objective=obj; row.message=output.message;
   row.boundary=abs(p-lb)<1e-5.*max(1,abs(lb)) | abs(p-ub)<1e-5.*max(1,abs(ub));
   row.boundary(~isfinite(ub))=abs(p(~isfinite(ub))-lb(~isfinite(ub)))<1e-5;
   row.jacobian_condition=cond(full(J));
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
  'fwhm_diagnostic',struct([]));
 if ~isempty(good)
  [~,ki]=min([candidates(good).objective]); best=good(ki); p=candidates(best).p;
  [~,bg,components]=prediction(p);
  pars=reshape(p(3:end),3,n).'; [~,order]=sort(pars(:,1)); pars=pars(order,:);
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
 if n==1, fits=f; else, fits(n)=f; end
end
 function [total,bg,comp]=prediction(p)
  bg=p(1)*(E/1000).^(-p(2)); comp=zeros(numel(E),(numel(p)-2)/3);
  for jj=1:size(comp,2)
   ii=3+(jj-1)*3;
   comp(:,jj)=model.model_fn(p(ii)*1000,p(ii+1)*1000,p(ii+2)*ampunit,E);
  end
  total=bg+sum(comp,2);
 end
end
