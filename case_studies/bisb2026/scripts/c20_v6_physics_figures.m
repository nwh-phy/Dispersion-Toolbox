function c20_v6_physics_figures(out)
a=load(fullfile(out,'appendix','products.mat')); products=a.products; cfg=a.cfg; rows={}; centroidrows={}; counterfactuals={};
for si=1:numel(products)
 if isempty(products{si}), continue; end
 product=products{si}; ds=product.details; map=product.map; bins=product.bins; sid=product.session;
 fig=figure('Visible','off','Position',[100 100 1150 750]);
 qi=abs(map.q)<=.015; ei=map.E>=200&map.E<=2100;
 image=imagesc(map.q(qi),map.E(ei)/1000,map.A1(ei,qi)); set(image,'AlphaData',isfinite(map.A1(ei,qi))); set(gca,'YDir','normal','Color',[.8 .8 .8]);
 colorbar; xlabel('Signed q (A^{-1})'); ylabel('Energy loss (eV)'); title([sid ' | A1 observed counts; no interpolation from selected members'],'Interpreter','none');
 hold on; yline(.3,'w--'); yline(1.8,'w--');
 for b=bins, if b.valid, plot([b.q_left b.q_right],[.3 .3],'r','LineWidth',3); end, end
 exportgraphics(fig,fullfile(out,'figures',[sid '_observed_map.png'])); close(fig);
 fig=figure('Visible','off','Position',[100 100 1150 750]);
 image=imagesc(map.q(qi),map.E(ei)/1000,map.A1(ei,qi)); set(image,'AlphaData',isfinite(map.A1(ei,qi))&map.A1(ei,qi)>0);
 set(gca,'YDir','normal','ColorScale','log','Color',[.8 .8 .8]); colorbar;
 xlabel('Signed q (A^{-1})'); ylabel('Energy (eV)'); title([sid ' | log color scale, original counts; nonpositive gray'],'Interpreter','none');
 exportgraphics(fig,fullfile(out,'figures',[sid '_observed_map_logcolor.png'])); close(fig);
 fig=figure('Visible','off','Position',[100 100 1300 850]); tiledlayout(2,ceil(numel(bins)/2),'TileSpacing','compact');
 for b=bins
  nexttile; plot(map.E/1000,b.mean,'k'); xlim([.3 1.8]); xlabel('eV'); ylabel('Counts (q mean)'); title(sprintf('q=%+.4f; N3',b.q_Ainv));
 end
 sgtitle([sid ' | observed spectra, independent ordinate scales'],'Interpreter','none');
 exportgraphics(fig,fullfile(out,'figures',[sid '_representative_spectra.png'])); close(fig);
 for di=1:numel(ds)
  d=ds{di}; f=d.fits(2); single=d.fits(1); q=d.unit.q_Ainv;
  if ~f.success, continue; end
  reference=[]; refq=sign(q)*cfg.frozen_reference_abs_q;
  for ri=1:numel(ds)
   z=ds{ri};
   if abs(z.unit.q_Ainv-refq)<1e-10&&strcmp(z.fits(2).peak_model,f.peak_model)&&z.fits(2).success
    reference=z.fits(2).components;
   end
  end
  [~,apex_indices]=max(f.components,[],1); apex_reversed=apex_indices(1)>apex_indices(2);
  if apex_reversed, reference=[]; end
  moments=qe_spectral_centroids(f.energy_meV,f.components,reference); assert(moments.identity_error<1e-8);
  observed_mu=trapz(f.energy_meV,f.energy_meV.*f.observed)/trapz(f.energy_meV,f.observed);
  [~,obspeak]=max(f.observed); [~,onepeak]=max(single.prediction-single.background);
  masklo=f.energy_meV<=900; maskhi=f.energy_meV>=900;
  low=trapz(f.energy_meV(masklo),f.observed(masklo)); high=trapz(f.energy_meV(maskhi),f.observed(maskhi));
  boundary=f.candidates(f.selected_start).boundary; background_boundary=any(boundary(1:2)); peak_boundary=any(boundary(3:end));
  for j=1:2
   [pk,apex]=max(f.components(:,j)); left=any(f.components(1:apex,j)<=pk/2); right=any(f.components(apex:end,j)<=pk/2);
   rows(end+1,:)={string(sid),q,string(f.peak_model),j,f.energy_meV(apex),f.parameters(j,1),f.parameters(j,2),f.parameters(j,3), ...
    moments.area(j),moments.fraction(j),moments.component_centroid(j),left&&right,background_boundary,peak_boundary,apex_reversed, ...
    string('conditional candidate; not a resolved mode or confidence interval')}; %#ok<AGROW>
  end
  centroidrows(end+1,:)={string(sid),q,string(f.peak_model),observed_mu,f.energy_meV(obspeak),single.energy_meV(onepeak), ...
   moments.combined_centroid,moments.frozen_centroid,moments.component_centroid(1),moments.component_centroid(2),moments.fraction(1), ...
   moments.identity_error,low,high,low/high,refq,moments.reference_applicable}; %#ok<AGROW>
  counterfactuals{end+1}=struct('session',sid,'q',q,'model',f.peak_model,'E',f.energy_meV,'moments',moments,'reference_q',refq); %#ok<AGROW>
 end
end
t=cell2table(rows,'VariableNames',{'session','q','model','component','apex_meV','native_E0_meV','native_width_meV','native_A','area_reference','fraction','centroid_meV', ...
 'half_height_bracketed','background_boundary','peak_boundary','apex_order_reversed','interpretation'});
c=cell2table(centroidrows,'VariableNames',{'session','q','model','observed_centroid','observed_grid_apex','single_component_grid_apex','two_component_centroid','frozen_shape_centroid', ...
 'mu1','mu2','f1','identity_error','observed_low_area','observed_high_area','observed_low_high_ratio','reference_q','reference_available'});
writetable(t,fullfile(out,'physics_candidates.csv')); writetable(c,fullfile(out,'spectral_weight_readout.csv')); save(fullfile(out,'frozen_shape_counterfactuals.mat'),'counterfactuals','-v7');
for sid=unique(t.session).'
 fig=figure('Visible','off','Position',[100 100 1200 780]); tiledlayout(2,2);
 for mi=1:2
  model=cfg.models{mi}; nexttile(mi); hold on;
  for j=1:2
   ix=t.session==sid&t.model==model&t.component==j;
   scatter(t.q(ix),t.apex_meV(ix)/1000,55,'DisplayName',sprintf('P%d conditional',j));
   bad=ix & (t.peak_boundary|t.apex_order_reversed); scatter(t.q(bad),t.apex_meV(bad)/1000,85,'x','HandleVisibility','off');
  end
  title([model ' | x = peak bound or reversed apex order'],'Interpreter','none'); xlabel('Signed q (A^{-1})'); ylabel('Component curve apex (eV)'); legend('show','Location','best'); ylim([.3 1.8]);
  nexttile(mi+2); hold on;
  for j=1:2, ix=t.session==sid&t.model==model&t.component==j; scatter(t.q(ix),t.fraction(ix),55,'DisplayName',sprintf('P%d',j)); end
  xlabel('Signed q (A^{-1})'); ylabel('W[0.3,1.8] fraction, excluding background'); ylim([0 1]);
 end
 sgtitle([char(sid) ' | templates are sensitivity, NOT confidence intervals'],'Interpreter','none');
 exportgraphics(fig,fullfile(out,'figures',[char(sid) '_candidate_scatter.png'])); close(fig);
 fig=figure('Visible','off','Position',[100 100 1150 800]); tiledlayout(2,2);
 for mi=1:2
  for signq=[-1 1]
   nexttile; ix=c.session==sid&c.model==cfg.models{mi}&sign(c.q)==signq;
   hold on; scatter(c.q(ix),c.observed_centroid(ix)/1000,45,'k','filled','DisplayName','Observed incl. background');
   scatter(c.q(ix),c.two_component_centroid(ix)/1000,55,'b','DisplayName','Actual two-component centroid');
   scatter(c.q(ix),c.frozen_shape_centroid(ix)/1000,60,'r','+','DisplayName','Frozen shapes, weights changed');
   xlabel('Signed q (A^{-1})'); ylabel('Window centroid (eV)'); title([cfg.models{mi} ' side ' num2str(signq)],'Interpreter','none'); legend('show','Location','best');
  end
 end
 sgtitle('Descriptive reweighting at fixed same-side |q|=0.0075; not causal partition');
 exportgraphics(fig,fullfile(out,'figures',[char(sid) '_centroid_reweighting.png'])); close(fig);
end
fig=figure('Visible','off','Position',[100 100 1200 750]); tiledlayout(2,2);
for mi=1:2
 for j=1:2
  nexttile; hold on;
  for si=1:2
   sid=string(products{si}.session); ix=t.session==sid&t.model==cfg.models{mi}&t.component==j;
   scatter(t.q(ix),t.apex_meV(ix)/1000,60,'DisplayName',sid);
  end
  xlabel('Signed q (A^{-1})'); ylabel(sprintf('P%d apex (eV)',j)); title(cfg.models{mi},'Interpreter','none'); legend('show','Interpreter','none');
 end
end
sgtitle('Frozen-method acquisition comparison; not blind or controlled thickness');
exportgraphics(fig,fullfile(out,'figures','repeat_comparison.png')); close(fig);
disp('PHYSICS_FIGURES_COMPLETE');
end
