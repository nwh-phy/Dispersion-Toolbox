function result = c20_v4_verify_packet(folder)
% Read only an extracted packet; independently reconstruct native model curves.
records={}; fitted=load(fullfile(folder,'590_PL2_10w','fit_details.mat'));
comparison=load(fullfile(folder,'590_PL2_10w','solver_background_comparison.mat'));
sets={fitted.details,comparison.independent,comparison.background}; nc=0; nm=0; npar=0;
for si=1:numel(sets)
 ds=sets{si};
 for j=1:numel(ds)
  d=ds{j};
  for f=d.fits
   nm=nm+1; npar=npar+f.n_components;
   peak=peak_models(f.peak_model); nc=nc+numel(f.candidates);
   for c=f.candidates
    if ~all(isfinite(c.p)), continue; end
    p=c.p; bg=p(1)*(f.energy_meV/1000).^(-p(2));
    if numel(p)>2+3*f.n_components, bg=bg+p(end); end
    curves=zeros(numel(f.energy_meV),f.n_components);
    native=reshape(p(3:2+3*f.n_components),3,[]).'.*[1000 1000 f.ampunit*f.scale];
    for k=1:f.n_components
     pp=native(k,:); curves(:,k)=peak.model_fn(pp(1),pp(2),pp(3),f.energy_meV);
    end
    pred=bg*f.scale+sum(curves,2); Q=sum((f.observed-pred).^2)/f.scale^2;
    assert(abs(Q-c.objective)<1e-7*max(1,Q),'Candidate objective mismatch');
   end
   if ~f.success, records(end+1,:)={si,string(d.key),f.n_components,false,NaN}; continue; end %#ok<AGROW>
   c=f.candidates(f.selected_start); p=c.p;
   native=reshape(p(3:2+3*f.n_components),3,[]).'.*[1000 1000 f.ampunit*f.scale];
   [~,order]=sort(native(:,1)); native=native(order,:);
   assert(isequal(order,f.component_order));
   inverse=zeros(1,f.n_components); inverse(order)=1:f.n_components;
   assert(isequal(inverse,f.raw_to_order));
   assert(max(abs(native-f.parameters),[],'all')<1e-7*max(1,max(abs(native),[],'all')));
   bg=p(1)*f.scale*(f.energy_meV/1000).^(-p(2));
   if numel(p)>2+3*f.n_components, bg=bg+p(end)*f.scale; end
   curves=zeros(size(f.components));
   for k=1:f.n_components
    pp=native(k,:); curves(:,k)=peak.model_fn(pp(1),pp(2),pp(3),f.energy_meV);
   end
   pred=bg+sum(curves,2); err=max(abs(pred-f.prediction))/max(1,max(abs(f.observed)));
   assert(err<1e-8); assert(max(abs(curves-f.components),[],'all')<1e-8*max(1,max(abs(f.observed))));
   assert(max(abs(f.observed-pred-f.residual))<1e-8*max(1,max(abs(f.observed))));
   assert(abs(sum(f.residual.^2)-f.sse)<1e-7*max(1,f.sse));
   if f.witness.available
    w=f.witness;
    assert(abs(sum((w.prediction-f.observed).^2)/f.scale^2-w.h0_objective)<1e-7*max(1,w.h0_objective));
   end
   records(end+1,:)={si,string(d.key),f.n_components,true,err}; %#ok<AGROW>
  end
 end
end
assert(nm==36&&npar==54&&nc==1008);
parent=readtable(fullfile(folder,'audit','boundary_by_parameter_native.csv'),TextType='string');
assert(height(parent)==16848&&nnz(parent.selected)==702);
assert(all(parent.flag_mismatch==0));
groups=findgroups(parent.key,parent.n_components,parent.start);
assert(max(groups)==2592);
for g=1:max(groups)
 b=parent(groups==g,:); n=b.n_components(1);
 raw=b.raw_value(3:2+3*n); raw=reshape(raw,3,n).';
 [~,order]=sort(raw(:,1)); inverse=zeros(1,n); inverse(order)=1:n;
 mapping=repelem(inverse,3).';
 assert(isequal(b.ordered_component(3:2+3*n),mapping));
 assert(max(abs(b.native_value-b.raw_value.*b.native_factor))<1e-8*max(1,max(abs(b.native_value))));
end
bins=load(fullfile(folder,'590_PL2_10w','centered_binned_spectra.mat'));
for b=bins.bins
 [ok,ix]=ismember(b.source_channel,bins.members); assert(all(ok));
 assert(max(abs(sum(bins.member_spectra(:,ix),2)-b.sum))<1e-8);
 assert(max(abs(b.sum/b.N-b.mean))<1e-8); assert(abs(b.q_right-b.q_left-b.q_width)<1e-12);
end
a=load(fullfile(folder,'590_PL2_10w','sequence_block_spectra.mat')); a=a.frames;
assert(size(a.block_sum_A0,2)==6&&size(a.block_sum_A0,3)==3);
for t=find(a.alignment.valid)
 idx=a.support+a.alignment.measured_offset_pixels(t);
 assert(isequaln(a.member_spectra_A1(:,t,:),a.member_spectra_A0(idx,t,:)));
end
for j=1:3
 [ok,ix]=ismember(a.n3_members{j},a.native_members); assert(all(ok));
 u0=mean(a.member_spectra_A0(a.support,:,ix),3); u1=mean(a.member_spectra_A1(:,:,ix),3);
 for b=1:6
  mask=a.block_id==b&a.alignment.valid;
  assert(max(abs(sum(u0(:,mask),2)-a.block_sum_A0(:,b,j)))<1e-8);
  assert(max(abs(sum(u1(:,mask),2)-a.block_sum_A1(:,b,j)))<1e-8);
  assert(max(abs(a.block_per_frame_mean_A0(:,b,j)*nnz(mask)-a.block_sum_A0(:,b,j)))<1e-8);
 end
end
result=struct('models_reconstructed',nm,'components_reconstructed',npar,'candidates_recomputed',nc, ...
 'parent_candidates_mapping_recomputed',2592, ...
 'bins_recomputed',numel(bins.bins),'real_block_spectra_recomputed',18, ...
 'A1_member_arrays_checked',true,'fit_checks',cell2table(records,'VariableNames',{'set','key','n_components','checked','relative_error'}));
disp('PACKET_ARRAY_READBACK_PASSED');
end
