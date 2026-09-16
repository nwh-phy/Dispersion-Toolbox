function stats = c20_v4_export_fits(details,dest)
% Reconstruct every selected curve and audit each candidate independently.
if ~isfolder(dest), mkdir(dest); end
br={}; cr={}; pr={}; mr={}; checks={}; wr={}; ar={}; maps=cell(size(details));
for di=1:numel(details)
 d=details{di}; u=d.unit; maps{di}=cell(1,numel(d.fits));
 for f=d.fits
  n=f.n_components; selected_ok=f.success&&isfinite(f.selected_start)&&f.selected_start>=1&&f.selected_start<=numel(f.candidates);
  cm=cell(1,numel(f.candidates)); any_selected=false;
  for k=1:numel(f.candidates)
   c=f.candidates(k); m=qe_component_mapping(c.p,f.lb,f.ub,n,f.scale,f.ampunit,f.peak_model); cm{k}=m;
   selected=selected_ok&&k==f.selected_start;
   hits=m.lower_hit|m.upper_hit;
   if selected, any_selected=any(hits); end
   cr(end+1,:)={string(d.key),u.bin_size_requested,u.source_q_count,n,k,selected,c.exitflag,c.objective,m.finite, ...
    any(hits),string(mat2str(m.order)),string(mat2str(m.raw_to_order)),string(mat2str(m.parameters,17)), ...
    m.equal_center,string(mat2str(c.p0,17)),string(mat2str(c.p,17))}; %#ok<AGROW>
   for j=1:numel(c.p)
    factor=m.native_scale(j);
    br(end+1,:)={string(d.key),u.bin_size_requested,u.source_q_count,n,u.q_Ainv,u.q_left,u.q_right, ...
     min(f.energy_meV),max(f.energy_meV),string(f.peak_model),k,selected,j,m.raw_component(j),m.ordered_component(j), ...
     m.names(j),m.units(j),c.p(j),f.lb(j),f.ub(j),factor,c.p(j)*factor,f.lb(j)*factor,f.ub(j)*factor, ...
     (c.p(j)-f.lb(j))*factor,(f.ub(j)-c.p(j))*factor,m.tol_lower(j),m.tol_upper(j), ...
     c.boundary(j),m.lower_hit(j),m.upper_hit(j),m.boundary_type(j),m.finite,c.exitflag, ...
     m.finite&&c.boundary(j)~=hits(j)}; %#ok<AGROW>
   end
  end
  maps{di}{n}=cm;
  mr(end+1,:)={string(d.key),u.bin_size_requested,n,u.q_Ainv,string(f.peak_model),f.success, ...
   f.selected_start,f.normalized_sse,any_selected,string(f.numerical_status)}; %#ok<AGROW>
  perr=NaN; cerr=NaN; rerr=NaN; qerr=NaN; parerr=NaN;
  if selected_ok
   c=f.candidates(f.selected_start); m=cm{f.selected_start};
   model=peak_models(f.peak_model); curves=zeros(size(f.components));
   for j=1:n
    p=m.parameters(j,:); curves(:,j)=model.model_fn(p(1),p(2),p(3),f.energy_meV);
   end
   bg=c.p(1)*f.scale*(f.energy_meV/1000).^(-c.p(2));
   if numel(c.p)>2+3*n, bg=bg+c.p(end)*f.scale; end
   pred=bg+sum(curves,2);
   parerr=max(abs(m.parameters-f.parameters),[],'all')/max(1,max(abs(f.parameters),[],'all'));
   perr=max(abs(pred-f.prediction))/max(1,max(abs(f.observed)));
   cerr=max(abs(curves-f.components),[],'all')/max(1,max(abs(f.observed)));
   rerr=max(abs(f.observed-pred-f.residual))/max(1,max(abs(f.observed)));
   qerr=abs(sum((f.observed-pred).^2)/f.scale^2-f.normalized_sse);
   assert(max([parerr perr cerr rerr qerr])<1e-7,'Saved fit contract mismatch');
  end
  checks(end+1,:)={string(d.key),n,selected_ok,parerr,perr,cerr,rerr,qerr}; %#ok<AGROW>
  for j=1:n
   area=NaN; area_ref=NaN; apex=NaN; width=NaN; bracketed=false; B=NaN; r=NaN; C=NaN;
   scope=qe_area_scopes(f.energy_meV,nan(size(f.energy_meV)));
   pars=f.parameters(j,:);
   if selected_ok
    curve=f.components(:,j); [pk,ip]=max(curve); apex=f.energy_meV(ip);
    left=find(curve(1:ip)<=pk/2,1,'last'); right=ip-1+find(curve(ip:end)<=pk/2,1,'first');
    bracketed=pk>0&&~isempty(left)&&~isempty(right)&&left<ip&&right>ip;
    if bracketed
     xl=interp1(curve(left:left+1),f.energy_meV(left:left+1),pk/2);
     xr=interp1(curve(right-1:right),f.energy_meV(right-1:right),pk/2); width=xr-xl;
    end
    scope=qe_area_scopes(f.energy_meV,curve);
    area=scope.area_fit_window; area_ref=scope.area_reference_window; p=f.candidates(f.selected_start).p;
    B=p(1)*f.scale; r=p(2); C=0; if numel(p)>2+3*n, C=p(end)*f.scale; end
   end
   pr(end+1,:)={string(d.key),n,j,u.q_Ainv,string(f.peak_model),pars(1),pars(2),pars(3), ...
    B,r,C,area,area_ref,apex,width,bracketed,string('not_assessed'),string('not_assessed'),string('not_assessed')}; %#ok<AGROW>
   ar(end+1,:)={string(mat2str(scope.fit_window_meV)),string(mat2str(scope.reference_window_meV)), ...
    string(mat2str(scope.reference_actual_support_meV)),scope.reference_window_fully_observed,string(scope.integration_method)}; %#ok<AGROW>
  end
  if isfield(f,'witness')&&f.witness.available
   w=f.witness; assert(abs(w.objective-w.h0_objective)<1e-7*max(1,w.h0_objective));
   wr(end+1,:)={string(d.key),w.h0_objective,f.normalized_sse,w.objective,w.best_feasible_objective, ...
    w.optimized_nesting_violation,string(w.candidate_type)}; %#ok<AGROW>
  end
 end
end
b=cell2table(br,'VariableNames',{'key','bin_N_requested','source_q_count','n_components','q_Ainv','q_left','q_right', ...
 'window_low','window_high','peak_model','start','selected','raw_parameter_index','raw_component','ordered_component', ...
 'parameter','native_unit','raw_value','raw_lb','raw_ub','native_factor','native_value','native_lb','native_ub', ...
 'native_distance_lower','native_distance_upper','legacy_tol_lower','legacy_tol_upper','raw_boundary','lower_hit','upper_hit', ...
 'boundary_type','finite_candidate','exitflag','flag_mismatch'});
c=cell2table(cr,'VariableNames',{'key','bin_N_requested','source_q_count','n_components','start','selected','exitflag','objective', ...
 'finite','boundary','order','raw_to_order','ordered_native_parameters','equal_center','p0','p'});
models=cell2table(mr,'VariableNames',{'key','bin_N_requested','n_components','q_Ainv','peak_model','success','selected_start','normalized_sse','boundary','numerical_status'});
pars=cell2table(pr,'VariableNames',{'key','n_components','ordered_component','q_Ainv','peak_model','E0_meV','native_width_meV','native_A', ...
 'background_B','background_r','background_C','area_fit_window','area_reference_window','apex_meV','zero_baseline_FWHM_meV','width_bracketed', ...
 'energy_identifiability','width_identifiability','area_identifiability'});
pars=[pars,cell2table(ar,'VariableNames',{'fit_window_meV','reference_window_meV','reference_actual_support_meV','reference_window_fully_observed','integration_method'})];
writetable(b,fullfile(dest,'boundary_by_parameter_native.csv'));
writetable(c,fullfile(dest,'candidate_solution_families.csv'));
writetable(models,fullfile(dest,'model_comparison.csv')); writetable(pars,fullfile(dest,'component_parameters.csv'));
writetable(cell2table(checks,'VariableNames',{'key','n_components','checked','parameter_error','prediction_error','component_error','residual_error','objective_error'}),fullfile(dest,'contract_checks.csv'));
types=unique(b.boundary_type(b.selected&(b.lower_hit|b.upper_hit))); tr=cell(0,3);
for t=types.'
 ix=b.selected & b.boundary_type==t;
 tr(end+1,:)={t,nnz(ix),numel(unique(b.key(ix)+"_n"+b.n_components(ix)))}; %#ok<AGROW>
end
writetable(cell2table(tr,'VariableNames',{'boundary_type','parameter_rows','distinct_models'}),fullfile(dest,'boundary_type_summary.csv'));
if ~isempty(wr), writetable(cell2table(wr,'VariableNames',{'key','Q1_optimized','Q2_optimized','Q2_witness','Q2_best_feasible','optimized_violation','witness_type'}),fullfile(dest,'nested_objective_checks.csv')); end
save(fullfile(dest,'component_mappings.mat'),'maps','-v7');
stats=struct('models',height(models),'components',height(pars),'candidates',height(c), ...
 'parameter_rows',height(b),'selected_parameter_rows',nnz(b.selected),'selected_boundary_models',nnz(models.boundary), ...
 'boundary_flag_mismatches',nnz(b.flag_mismatch),'failed_models',nnz(~models.success));
c20_v4_io('json',fullfile(dest,'inventory.json'),stats);
end
