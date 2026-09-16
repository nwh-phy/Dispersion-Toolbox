function result = c20_v5_verify(folder)
% Independent readback from MAT observations and native model definitions.
a=load(fullfile(folder,'A0_A1','fit_details.mat')); A1=a.A1;
a=load(fullfile(folder,'A0_A1','A0_existing_mode2.mat')); A0=a.A0;
s=load(fullfile(folder,'sequence_models','sequence_inputs.mat')); frames=s.frames;
nm=0; nc=0;
for di=1:numel(A1)
 d=A1{di}; j=ceil(di/2);
 for f=d.fits
  nm=nm+1; nc=nc+numel(f.candidates); model=peak_models(f.peak_model);
  mask=frames.E>=min(f.energy_meV)&frames.E<=max(f.energy_meV);
  y1=sum(frames.sequence_bin_A1(mask,:,j),2); y0=sum(frames.sequence_bin_A0(mask,:,j),2);
  assert(max(abs(y1-f.observed))<1e-8); assert(max(abs(y0-A0{di}.fits(f.n_components).observed))<1e-8);
  for c=f.candidates
   if ~all(isfinite(c.p)), continue; end
   p=c.p; pred=p(1)*f.scale*(f.energy_meV/1000).^(-p(2));
   for k=1:f.n_components
    ii=3+(k-1)*3; pred=pred+model.model_fn(p(ii)*1000,p(ii+1)*1000,p(ii+2)*f.ampunit*f.scale,f.energy_meV);
   end
   assert(abs(sum((pred-f.observed).^2)/f.scale^2-c.objective)<1e-7*max(1,c.objective));
  end
  if f.success
   C=zeros(size(f.components));
   for k=1:f.n_components, pp=f.parameters(k,:); C(:,k)=model.model_fn(pp(1),pp(2),pp(3),f.energy_meV); end
   assert(max(abs(C-f.components),[],'all')<1e-8*max(1,max(abs(f.observed))));
   assert(max(abs(f.background+sum(C,2)-f.prediction))<1e-8*max(1,max(abs(f.observed))));
   assert(max(abs(f.observed-f.prediction-f.residual))<1e-8);
  end
 end
end
assert(nm==12&&nc==360);
m=load(fullfile(folder,'member_models','all_member_fits.mat')); member_nc=0;
for di=1:numel(m.allmembers)
 d=m.allmembers{di}; region=d.region; idx=(region-1)*5+(1:5);
 for f=d.fits
  assert(f.observation_count==numel(f.observed)); Q=numel(f.q); member_nc=member_nc+numel(f.candidates); model=peak_models(f.peak_model);
  emask=frames.native_E>=min(f.E)&frames.native_E<=max(f.E);
  observed=squeeze(sum(frames.member_spectra_A0(emask,:,idx),2)); assert(max(abs(observed-f.observed),[],'all')<1e-8);
  beta=(f.q-min(f.q))/(max(f.q)-min(f.q));
  for c=f.candidates
   if ~all(isfinite(c.p)), continue; end
   p=c.p; pred=(f.E/1000).^(-p(1))*p(2:Q+1)*f.scale;
   for k=1:f.n_components
    b=1+Q+(k-1)*(3+Q); centers=1000*((1-beta)*p(b+1)+beta*p(b+2));
    for q=1:Q
     pred(:,q)=pred(:,q)+model.model_fn(centers(q),p(b+3)*1000,p(b+3+q)*f.ampunit*f.scale,f.E);
    end
   end
   assert(abs(sum((pred-f.observed).^2,'all')/f.scale^2-c.objective)<1e-7*max(1,c.objective));
  end
  if f.success
   C=zeros(size(f.components));
   for k=1:f.n_components
    for q=1:Q, C(:,q,k)=model.model_fn(f.native_centers(q,k),f.native_widths(k),f.native_A(q,k),f.E); end
   end
   assert(max(abs(C-f.components),[],'all')<1e-8*max(1,max(abs(f.observed),[],'all')));
   assert(max(abs(sum(C,3)+f.background-f.prediction),[],'all')<1e-8*max(1,max(abs(f.observed),[],'all')));
  end
  if f.witness.available
   assert(max(abs(f.witness.prediction-d.fits(1).prediction),[],'all')<1e-8);
  end
 end
 b=load(fullfile(folder,'member_models',d.key,'derived_bin_predictions.mat'));
 for n=1:2
  for v=b.predictions{n}
   [ok,ix]=ismember(v.members,d.members); assert(all(ok));
   assert(max(abs(v.predicted_mean-mean(d.fits(n).prediction(:,ix),2)))<1e-8);
  end
 end
end
areas=load(fullfile(folder,'audit','area_arrays.mat')); count=0;
for cell=areas.area_records
 rec=cell{1}; ref=rec.E>=300&rec.E<=1800;
 assert(max(abs(trapz(rec.E(ref),rec.components(ref,:),1)-rec.areas.area_reference_window))<1e-8);
 assert(max(abs(trapz(rec.E,rec.components,1)-rec.areas.area_fit_window))<1e-8); count=count+rec.n;
end
assert(count==162&&member_nc==144);
result=struct('A1_models',nm,'A1_candidates',nc,'member_models',12,'member_candidates',member_nc, ...
 'area_rows_recomputed',count,'member_arrays_and_derived_bins',true,'same_support_A0_A1',true);
disp('V5_ARRAY_READBACK_PASSED');
end
