function c20_v5_qsequence(out)
parent=jsondecode(fileread(fullfile(out,'parent_packet_reference.json'))); parent=parent.parent;
p=jsondecode(fileread(fullfile(parent,'provenance','parents.json')));
manifest=jsondecode(fileread(fullfile(p.v2,'input_manifest.resolved.yaml'))); rec=manifest.sessions(strcmp({manifest.sessions.session_id},'590_PL2_10w'));
rawfile=rec.files(find(endsWith(string({rec.files.path}),'.npy'),1));
assert(strcmpi(c20_v4_io('hash',rawfile.path),rawfile.sha256));
a=load(fullfile(parent,'590_PL2_10w','sequence_block_spectra.mat')); frames=a.frames;
raw=read_npy(rawfile.path); assert(isequal(size(raw),[300 512 1028]));
center=rec.actual.q_zero_native_channel; members=center+(-10:10); assert(~any(ismember(members,rec.actual.invalid_native_q)));
idx=find(frames.native_E>=-100&frames.native_E<=100); profiles_A0=zeros(21,300); profiles_A1=profiles_A0;
for t=1:300
 profiles_A0(:,t)=sum(double(squeeze(raw(t,members,idx))),2);
 profiles_A1(:,t)=sum(double(squeeze(raw(t,members,idx+frames.alignment.measured_offset_pixels(t)))),2);
end
clear raw
q=((members-center)*rec.actual.dq_Ainv).'; rows={}; centers=zeros(300,6);
for t=1:300
 [~,ip]=max(profiles_A1(:,t));
 for band=1:3
  half=[2 4 8]; ix=abs(members-center)<=half(band); x=members(ix).'; y=profiles_A1(ix,t);
  centroid=sum(x.*y)/sum(y); baseline=mean(y([1 end])); yc=max(0,y-baseline);
  corrected=sum(x.*yc)/sum(yc); centers(t,2*band-1:2*band)=[centroid corrected];
  rows(end+1,:)={t,half(band),members(ip),centroid,corrected,baseline}; %#ok<AGROW>
 end
end
dest=fullfile(out,'q_center_diagnostics');
writetable(cell2table(rows,'VariableNames',{'sequence','half_band_channels','argmax_channel','raw_centroid_channel','edge_subtracted_centroid_channel','edge_background'}),fullfile(dest,'centroid_sensitivity.csv'));
save(fullfile(dest,'elastic_profiles.mat'),'profiles_A0','profiles_A1','members','q','centers','idx','-v7');
figure1=figure('Visible','off','Position',[100 100 1100 800]); tiledlayout(2,1);
nexttile; plot(1:300,centers); hold on; [~,arg]=max(profiles_A1,[],1); stairs(1:300,members(arg),'k:');
ylabel('Native channel'); xlabel('Sequence frame'); title('Elastic q profile: centroid band sensitivity; no q realignment');
nexttile; hold on;
for b=1:6
 y=mean(profiles_A1(:,frames.block_id==b),2); plot(members,y/sum(y),'DisplayName',sprintf('Block %d',b));
end
legend('show'); xlabel('Native channel'); ylabel('Display-only normalized elastic profile');
exportgraphics(figure1,fullfile(dest,'q_profiles.png')); close(figure1);
c20_v4_io('json',fullfile(dest,'status.json'),struct('q_center_stability','unresolved: centroid/argmax changes may include beam profile deformation', ...
 'q_realignment','not_performed','reference_band_channels',members,'energy_window',[-100 100],'raw_sha256',rawfile.sha256));
rows={}; block_arrays=struct('E',frames.E,'A0',frames.block_per_frame_mean_A0,'A1',frames.block_per_frame_mean_A1,'count',frames.block_count,'q',frames.targets);
E=frames.E; mask=E>=300&E<=1800; lo=E>=300&E<=900; hi=E>=900&E<=1800;
for j=1:3
 reference=mean(frames.block_per_frame_mean_A1(mask,:,j),2);
 fig=figure('Visible','off','Position',[100 100 1100 800]); tiledlayout(2,1); nexttile; hold on;
 for b=1:6
  yy=frames.block_per_frame_mean_A1(:,b,j); y=yy(mask); coefficient=(reference.'*y)/(reference.'*reference);
  residue=y-coefficient*reference; area=trapz(E(mask),y); energy_center=trapz(E(mask),E(mask).*y)/area;
  rows(end+1,:)={j,b,frames.targets(j),area,energy_center,trapz(E(lo),yy(lo)),trapz(E(hi),yy(hi)),coefficient,sqrt(mean(residue.^2))}; %#ok<AGROW>
  plot(E(mask),y,'DisplayName',sprintf('Block %d',b));
 end
 legend('show'); ylabel('Counts / frame'); title(sprintf('R%d A1 same-region continuous sequence; independence not assumed',j));
 nexttile; hold on;
 for b=1:6
  y=frames.block_per_frame_mean_A1(mask,b,j); coefficient=(reference.'*y)/(reference.'*reference); plot(E(mask),y-coefficient*reference);
 end
 xlabel('Energy (meV)'); ylabel('Residual from scale-only block model');
 exportgraphics(fig,fullfile(out,'sequence_models',sprintf('blocks_R%d.png',j))); close(fig);
end
writetable(cell2table(rows,'VariableNames',{'region','block','q','area_300_1800_per_frame','spectral_centroid_meV','low_300_900_area','high_900_1800_area','scale_only_coefficient','scale_only_residual_RMS'}),fullfile(out,'sequence_models','block_diagnostics.csv'));
save(fullfile(out,'sequence_models','block_diagnostics.mat'),'block_arrays','-v7');
c20_v4_io('json',fullfile(out,'sequence_models','status.json'),struct('same_region','user_confirmed','continuous','user_confirmed','scan_move_adjust','user_confirmed_absent', ...
 'stationarity','not_established','block_fits','not_run: descriptive spectra and scale-only comparison only','independent_blocks',false));
assert(strcmpi(c20_v4_io('hash',rawfile.path),rawfile.sha256)); disp('Q_AND_SEQUENCE_DIAGNOSTICS_COMPLETE');
end
