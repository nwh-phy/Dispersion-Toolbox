function c20_v4_profiles(out)
% Descriptive objective profiles only; no likelihood confidence thresholds.
a=load(fullfile(out,'590_PL2_10w','solver_background_comparison.mat'));
dest=fullfile(out,'validation','profiles'); assert(~isfolder(dest)); mkdir(dest); rows={}; profiles={};
opt=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200, ...
 'MaxFunctionEvaluations',4000,'ConstraintTolerance',1e-7,'OptimalityTolerance',1e-7);
for j=1:3
 d=a.independent{2*j-1}; f=d.fits(2); assert(f.success);
 best=f.candidates(f.selected_start).p; m=qe_component_mapping(best,f.lb,f.ub,2,f.scale,f.ampunit,f.peak_model);
 ps=reshape(best(3:8),3,2).'; ps=ps(m.order,:); best(3:8)=reshape(ps.',1,6);
 [~,~,comp]=qe_component_prediction(f.energy_meV,best,f.peak_model,2,f.ampunit);
 areas=trapz(f.energy_meV,comp); ratio=areas(2)/sum(areas);
 grids={unique([0 (best(6)-best(3))*1000 900]),unique([4 best(4)*1000 5000]),unique([.05 ratio .95])};
 names={'separation_meV','P1_width_meV','P2_reference_area_fraction'};
 for variable=1:3
  for value=grids{variable}
   lb=f.lb; ub=f.ub; Aeq=[]; beq=[]; nonlinear=[];
   if variable==1, Aeq=zeros(1,8); Aeq(6)=1; Aeq(3)=-1; beq=value/1000; end
   if variable==2, lb(4)=value/1000; ub(4)=value/1000; end
   if variable==3, nonlinear=@(p)ratio_constraint(p,value); end
   trials=struct([]);
   for k=1:3
    p0=best;
    if k==2, p0([4 7])=min(5,p0([4 7])*1.4); end
    if k==3, p0([5 8])=p0([5 8])*.5; end
    p0=max(lb,min(ub,p0));
    order_constraint=zeros(1,8); order_constraint(3)=1; order_constraint(6)=-1;
    [p,Q,flag,output]=fmincon(@objective,p0,order_constraint,0,Aeq,beq,lb,ub,nonlinear,opt);
    row=struct('p',p,'Q',Q,'exitflag',flag,'output',output);
    if isempty(trials), trials=row; else, trials(end+1)=row; end
   end
   valid=[trials.exitflag]>0; Q=NaN;
   if any(valid), Q=min([trials(valid).Q]); end
   rows(end+1,:)={j,string(names{variable}),value,Q,nnz(valid),string('objective_not_confidence_interval')}; %#ok<AGROW>
   profiles{end+1}=struct('q',d.unit.q_Ainv,'variable',names{variable},'value',value,'trials',trials,'fit',f); %#ok<AGROW>
  end
 end
 disp(sprintf('PROFILE R%d complete',j));
end
writetable(cell2table(rows,'VariableNames',{'target','parameter','fixed_value','best_converged_objective','converged_starts','interpretation'}),fullfile(dest,'profile_objective.csv'));
save(fullfile(dest,'profile_arrays.mat'),'profiles','-v7');
 function q=objective(p)
  y=qe_component_prediction(f.energy_meV,p,f.peak_model,2,f.ampunit); q=sum((y-f.observed/f.scale).^2);
 end
 function [c,ceq]=ratio_constraint(p,value)
  [~,~,curves]=qe_component_prediction(f.energy_meV,p,f.peak_model,2,f.ampunit);
  ar=trapz(f.energy_meV,curves); c=[]; ceq=ar(2)/max(sum(ar),eps)-value;
 end
end
