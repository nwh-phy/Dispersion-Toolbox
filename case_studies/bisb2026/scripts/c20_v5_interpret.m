function c20_v5_interpret(out)
a=load(fullfile(out,'A0_A1','A0_existing_mode2.mat')); A0=a.A0; a=load(fullfile(out,'A0_A1','fit_details.mat')); A1=a.A1;
rows={};
for di=1:numel(A1)
 for n=1:2
  f0=A0{di}.fits(n); f1=A1{di}.fits(n); area0=qe_area_scopes(f0.energy_meV,f0.components); area1=qe_area_scopes(f1.energy_meV,f1.components);
  for j=1:n
   [~,k0]=max(f0.components(:,j)); [~,k1]=max(f1.components(:,j));
   rows(end+1,:)={ceil(di/2),string(f1.peak_model),n,j,f0.parameters(j,1),f1.parameters(j,1),f1.parameters(j,1)-f0.parameters(j,1), ...
    f0.energy_meV(k0),f1.energy_meV(k1),f0.parameters(j,2),f1.parameters(j,2),area0.area_reference_window(j),area1.area_reference_window(j), ...
    area0.area_reference_window(j)/sum(area0.area_reference_window),area1.area_reference_window(j)/sum(area1.area_reference_window)}; %#ok<AGROW>
  end
 end
end
writetable(cell2table(rows,'VariableNames',{'region','model','n','component','E0_A0','E0_A1','delta_E0','apex_A0','apex_A1','width_A0','width_A1','area_A0','area_A1','fraction_A0','fraction_A1'}),fullfile(out,'A0_A1','parameter_changes.csv'));
a=load(fullfile(out,'member_models','all_member_fits.mat')); rows={};
for cell=a.allmembers
 d=cell{1}; m1=d.fits(1); m2=d.fits(2);
 rows(end+1,:)={d.region,string(m1.peak_model),m1.objective,m2.objective,(m1.objective-m2.objective)/m1.objective, ...
  m1.slope_meV_A(1),string(mat2str(m2.slope_meV_A)),m1.native_widths(1),string(mat2str(m2.native_widths)), ...
  sqrt(mean(movmean(m1.residual,9,1).^2,'all')),sqrt(mean(movmean(m2.residual,9,1).^2,'all')), ...
  m2.witness.optimized_violation}; %#ok<AGROW>
end
writetable(cell2table(rows,'VariableNames',{'region','model','Q_M1','Q_M2','relative_gain_descriptive','M1_slope_meV_A','M2_slopes_meV_A','M1_width_meV','M2_widths_meV','M1_smooth_residual_RMS','M2_smooth_residual_RMS','nesting_violation'}),fullfile(out,'member_models','competition_summary.csv'));
end
