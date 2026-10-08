function frames = c20_v4_frames(rec,qe,bins,dest)
% Read native 590 once, retain only small representative member arrays.
files=rec.files; rawfile=files(find(endsWith(string({files.path}),'.npy'),1));
jsonfile=files(find(endsWith(string({files.path}),'.json'),1));
assert(strcmpi(c20_v4_io('hash',rawfile.path),rawfile.sha256),'Raw hash mismatch');
assert(strcmpi(c20_v4_io('hash',jsonfile.path),jsonfile.sha256),'JSON hash mismatch');
meta=jsondecode(fileread(jsonfile.path)); raw=read_npy(rawfile.path);
assert(meta.is_sequence&&meta.datum_dimension_count==2&&meta.collection_dimension_count==0);
assert(size(raw,2)==numel(qe.q_Ainv)&&size(raw,3)==numel(qe.energy_meV),'Explicit sequence-q-E axes mismatch');
T=size(raw,1); nq=size(raw,2); ne=size(raw,3); invalid=false(1,nq); a0=zeros(ne,nq);
badcounts=zeros(1,nq); observed_max=zeros(1,nq);
for q=1:nq
 x=double(squeeze(raw(:,q,:))); bad=~isfinite(x)|x>=double(intmax('uint32'));
 invalid(q)=any(bad,'all'); badcounts(q)=nnz(bad); observed_max(q)=max(x,[],'all');
 if ~invalid(q), a0(:,q)=sum(x,1).'; else, a0(:,q)=NaN; end
end
assert(isequal(find(invalid),rec.actual.invalid_native_q(:).'),'Parent mask mismatch');
identity_error=max(abs(a0(:,~invalid)-qe.intensity(:,~invalid)),[],'all');
assert(identity_error<1e-8*max(a0(:,~invalid),[],'all'),'A0 parent sum mismatch');
members=unique([bins.source_channel]); X=permute(double(raw(:,members,:)),[3 1 2]);
refq=find(abs(qe.q_Ainv)<=.001 & ~invalid);
reference=squeeze(sum(double(raw(:,refq,:)),2)).';
integrated=zeros(T,1); qcenter=nan(T,1);
for t=1:T
 x=double(squeeze(raw(t,~invalid,:))); integrated(t)=sum(x,'all');
 z=sum(x(:,abs(qe.energy_meV)<=100),2); validq=find(~invalid); [~,ii]=max(z); qcenter(t)=validq(ii);
end
clear raw
alignment=qe_zlp_integer_align(qe.energy_meV,reference,X,[-100 100]);
ar=qe_zlp_integer_align(qe.energy_meV,reference,reference,[-100 100]);
support=alignment.support; n3=bins([bins.N]==3); block_id=min(6,1+floor((0:T-1)*6/T));
block_sum=zeros(numel(support),6,numel(n3)); block_sum_a1=block_sum;
sequence_a0=zeros(numel(support),T,numel(n3)); sequence_a1=sequence_a0;
block_count=zeros(6,1);
for j=1:numel(n3)
 [ok,ix]=ismember(n3(j).source_channel,members); assert(all(ok));
 sequence_a0(:,:,j)=mean(X(support,:,ix),3);
 sequence_a1(:,:,j)=mean(alignment.aligned(:,:,ix),3);
 for b=1:6
  mask=block_id==b & alignment.valid; block_count(b)=nnz(mask);
  block_sum(:,b,j)=sum(sequence_a0(:,mask,j),2);
  block_sum_a1(:,b,j)=sum(sequence_a1(:,mask,j),2);
 end
end
frames=struct('E',alignment.E,'native_E',qe.energy_meV,'native_members',members, ...
 'q_members',qe.q_Ainv(members),'member_spectra_A0',X,'member_spectra_A1',alignment.aligned, ...
 'support',support,'reference_profiles_A0',reference,'reference_profiles_A1',ar.aligned, ...
 'reference_q_channels',qe.source_channel(refq),'invalid_q',invalid,'A0_identity_error',identity_error, ...
 'alignment',rmfield(alignment,'aligned'),'block_id',block_id,'block_count',block_count, ...
 'block_sum_A0',block_sum,'block_sum_A1',block_sum_a1, ...
 'block_per_frame_mean_A0',block_sum./reshape(block_count,1,6,1), ...
 'block_per_frame_mean_A1',block_sum_a1./reshape(block_count,1,6,1), ...
 'sequence_bin_A0',sequence_a0,'sequence_bin_A1',sequence_a1, ...
 'targets',[n3.q_Ainv],'n3_members',{ {n3.source_channel} }, ...
 'integrated_signal',integrated,'q_center_channel',qcenter);
save(fullfile(dest,'sequence_block_spectra.mat'),'frames','-v7');
h=meta.metadata.hardware_source;
semantics=struct('native_shape',[T nq ne],'axis_order','sequence,q,energy', ...
 'is_sequence',meta.is_sequence,'nimages',h.detector_configuration.nimages, ...
 'metadata_exposure_seconds',h.exposure,'frame_time_seconds',h.detector_configuration.frame_time, ...
 'sequence_semantics','camera sequence supported; fixed position and per-frame timestamps not verified', ...
 'same_location_and_stationarity','unknown','noise','not calibrated; detector corrections and sequence trends', ...
 'P5_role','engineering simulations only; no experimental false-positive rate', ...
 'raw_sha256',rawfile.sha256,'json_sha256',jsonfile.sha256);
c20_v4_io('json',fullfile(dest,'frame_semantics.json'),semantics);
writetable(table((1:T)',alignment.peak_pixel',alignment.valid',alignment.measured_offset_pixels',alignment.correction_pixels', ...
 integrated,qcenter,'VariableNames',{'sequence_index','ZLP_peak_pixel','ZLP_valid','measured_offset_pixels','correction_pixels','integrated_signal','q_center_channel'}),fullfile(dest,'frame_qc.csv'));
writetable(table((1:nq)',invalid',badcounts',observed_max','VariableNames',{'native_q','invalid','invalid_samples','max_native_value'}),fullfile(dest,'detector_mask.csv'));
f=figure('Visible','off','Position',[100 100 1000 750]); tiledlayout(3,1);
nexttile; plot(1:T,alignment.measured_offset_pixels*median(diff(qe.energy_meV))); ylabel('ZLP offset (meV)');
nexttile; plot(1:T,integrated); ylabel('All-q signal (counts)');
nexttile; plot(1:T,qcenter); ylabel('Central q channel'); xlabel('Sequence index (not independent repeats)');
exportgraphics(f,fullfile(dest,'figures','sequence_QC.png')); close(f);
for j=1:numel(n3)
 f=figure('Visible','off','Position',[100 100 1100 800]); tiledlayout(2,1);
 nexttile; plot(frames.E,frames.block_per_frame_mean_A0(:,:,j)); xlim([300 1800]); ylabel('Counts per sequence frame'); title(sprintf('A0 six actual blocks; q=%+.4f',n3(j).q_Ainv)); legend(compose('Block %d',1:6));
 nexttile; plot(frames.E,sum(sequence_a0(:,:,j),2),'k'); hold on; plot(frames.E,sum(sequence_a1(:,:,j),2),'r'); xlim([300 1800]); legend('A0','A1 integer common shift'); ylabel('Sequence-summed q mean'); xlabel('Energy (meV)');
 exportgraphics(f,fullfile(dest,'figures',sprintf('blocks_A0_A1_R%d.png',j))); close(f);
end
f=figure('Visible','off','Position',[100 100 1000 550]);
plot(qe.energy_meV,sum(reference,2),'k'); hold on; plot(frames.E,sum(ar.aligned,2),'r'); xlim([-100 100]); legend('A0 elastic reference','A1 elastic reference'); xlabel('Energy (meV)'); ylabel('Counts');
exportgraphics(f,fullfile(dest,'figures','ZLP_A0_A1.png')); close(f);
assert(strcmpi(c20_v4_io('hash',rawfile.path),rawfile.sha256));
end
