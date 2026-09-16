function c20_v5_export_members(d,dest)
if ~isfolder(dest), mkdir(dest); end
save(fullfile(dest,'member_fit.mat'),'d','-v7');
params={}; candidates={}; metrics={}; predictions=cell(1,2);
for f=d.fits
 n=f.n_components; predictions{n}=struct([]);
 for k=1:numel(f.candidates)
  c=f.candidates(k); opt=NaN; if isfield(c.output,'firstorderopt'), opt=c.output.firstorderopt; end
  candidates(end+1,:)={n,k,c.exitflag,c.objective,opt,string(mat2str(c.p0,17)),string(mat2str(c.p,17)),string(mat2str(c.order)),string(mat2str(c.boundary))}; %#ok<AGROW>
 end
 for j=1:n
  for q=1:numel(f.q)
   curve=f.components(:,q,j); a=qe_area_scopes(f.E,curve); [~,ip]=max(curve);
   params(end+1,:)={n,j,d.members(q),f.q(q),f.native_centers(q,j),f.native_widths(j),f.native_A(q,j), ...
    f.slope_meV_A(j),a.area_fit_window,a.area_reference_window,f.E(ip)}; %#ok<AGROW>
  end
 end
 for N=[1 3 5]
  center=ceil(numel(f.q)/2); ix=center-(N-1)/2:center+(N-1)/2;
  p=struct('N',N,'members',d.members(ix),'observed_mean',mean(f.observed(:,ix),2), ...
   'predicted_mean',mean(f.prediction(:,ix),2),'predicted_sum',sum(f.prediction(:,ix),2), ...
   'component_mean',squeeze(mean(f.components(:,ix,:),2)),'background_mean',mean(f.background(:,ix),2));
  if isempty(predictions{n}), predictions{n}=p; else, predictions{n}(end+1)=p; end
 end
 for q=1:numel(f.q)
  r=f.residual(:,q); y=f.observed(:,q);
  metrics(end+1,:)={n,d.members(q),f.q(q),sqrt(mean(r.^2)),sqrt(mean(movmean(r,9).^2)), ...
   sqrt(mean(y.^2)),f.success,f.objective}; %#ok<AGROW>
 end
end
writetable(cell2table(params,'VariableNames',{'n_components','center_ordered_component','native_member','q','E0_meV','width_meV','native_A','slope_meV_A', ...
 'area_fit_window','area_reference_window','apex_meV'}),fullfile(dest,'parameters.csv'));
writetable(cell2table(candidates,'VariableNames',{'n_components','start','exitflag','objective','firstorderopt','p0','p','order','boundary'}),fullfile(dest,'candidates.csv'));
writetable(cell2table(metrics,'VariableNames',{'n_components','member','q','residual_RMS','moving9_residual_RMS','observed_RMS','success','region_objective'}),fullfile(dest,'member_residual_diagnostics.csv'));
save(fullfile(dest,'derived_bin_predictions.mat'),'predictions','-v7');
fig=figure('Visible','off','Position',[100 100 1500 780]); tiledlayout(3,5,'TileSpacing','compact');
for q=1:5
 nexttile(q); plot(d.fits(1).E,d.fits(1).observed(:,q),'Color',[.7 .7 .7]); hold on;
 plot(d.fits(1).E,d.fits(1).prediction(:,q),'b'); plot(d.fits(2).E,d.fits(2).prediction(:,q),'r');
 title(sprintf('q=%+.4f; col %d',d.fits(1).q(q),d.members(q))); if q==1, legend('obs','M1','M2'); end
 nexttile(5+q); plot(d.fits(1).E,d.fits(1).residual(:,q),'b'); yline(0,'k:'); if q==1, ylabel('M1 residual'); end
 nexttile(10+q); plot(d.fits(2).E,d.fits(2).residual(:,q),'r'); yline(0,'k:'); xlabel('meV'); if q==1, ylabel('M2 residual'); end
end
sgtitle([d.key ' | observed five members only; local linear centers, shared widths'],'Interpreter','none');
exportgraphics(fig,fullfile(dest,'member_residuals.png')); close(fig);
end
