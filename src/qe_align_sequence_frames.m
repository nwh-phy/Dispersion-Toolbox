function A = qe_align_sequence_frames(source, zero, dq, opts)
%QE_ALIGN_SEQUENCE_FRAMES  Per-frame ZLP alignment and beam drift of a raw q-EELS sequence.
%   A = qe_align_sequence_frames(npy_path, zero, dq) reads a Nion Swift
%   sequence (frames x q x E .npy, calibration in the .json next to it),
%   builds the energy axis from the JSON and q = (channel - zero) * dq, and
%   aligns every frame by one integer energy shift (qe_zlp_integer_align)
%   measured on the channels |q| <= reference_q within alignment_window around
%   the ZLP found in raw_zlp_search. Only channels with |q| <= q_select are kept.
%   Then, per frame, the beam position along q is the centroid of the elastic
%   line (|E| <= q_centroid_energy, |q| <= q_centroid_range); drift_Ainv is its
%   q_median_frames moving median about the median (zero if do_q_drift = false).
%   beam is the per-frame reference count within |E| <= 100 meV, sigma_probe the
%   angular spread of the elastic line from the central five channels.
%
%   This is the raw-frame preparation of b1_low_branch_matlab, kept in one
%   place for every analysis of raw sequences. Bin the aligned frames with
%   qe_transport_frames, which follows drift_Ainv.
%
%   Fields of A:
%     source, raw_size        input file and its size (frames x q x E)
%     E                       energy (meV) relative to the reference ZLP pixel
%     Q, selected             q (1/A) and native channel of every kept channel
%     F                       aligned counts, E x frame x channel (double)
%     alignment               qe_zlp_integer_align output without 'aligned'
%     beam                    1 x T reference counts
%     centroid_Ainv           T x 1 elastic centroid (unsmoothed)
%     drift_Ainv              1 x T smoothed drift used for binning
%     sigma_probe             1/A
%     energy_origin_meV       raw energy of the reference ZLP pixel
%
%   See also: qe_transport_frames, qe_zlp_integer_align, read_npy
arguments
    source (1,:) char {mustBeFile}
    zero (1,1) double
    dq (1,1) double {mustBePositive}
    opts.q_select (1,1) double = 0.022
    opts.reference_q (1,1) double = 0.001
    opts.raw_zlp_search (1,2) double = [-350 100]
    opts.alignment_window (1,2) double = [-100 100]
    opts.q_centroid_range (1,1) double = 0.002
    opts.q_centroid_energy (1,1) double = 16
    opts.q_median_frames (1,1) double = 11
    opts.do_q_drift (1,1) logical = true
end
raw=read_npy(source);
assert(ndims(raw)==3,'qe_align_sequence_frames:Axes','Expected a frames x q x E sequence');
raw_size=size(raw);
meta=jsondecode(fileread(strrep(source,'.npy','.json')));ax=meta.spatial_calibrations;
eidx=find(contains(string({ax.units}),'eV'),1,'last');
Eraw=(ax(eidx).offset+(0:size(raw,3)-1)*ax(eidx).scale)'*1000;
q=((1:size(raw,2))-zero)*dq;refq=abs(q)<=opts.reference_q+1e-12;
ref=squeeze(sum(double(raw(:,refq,:)),2)).';
search=find(Eraw>=opts.raw_zlp_search(1)&Eraw<=opts.raw_zlp_search(2));
[~,p]=max(sum(ref(search,:),2));origin=Eraw(search(p));
selected=find(abs(q)<=opts.q_select+1e-12);
X=permute(raw(:,selected,:),[3 1 2]);clear raw
assert(all(isfinite(X(:)))&&all(X(:)<double(intmax('uint32'))),'qe_align_sequence_frames:Samples', ...
    'Invalid detector samples');
al=qe_zlp_integer_align(Eraw-origin,ref,X,opts.alignment_window);clear X
assert(all(al.valid),'qe_align_sequence_frames:Frames','Invalid ZLP frames require review');
E=Eraw(al.support)-Eraw(al.reference_pixel);T=numel(al.valid);Q=q(selected);
F=al.aligned;al=rmfield(al,'aligned');
beam=zeros(1,T);
for t=1:T,beam(t)=sum(ref(al.support(abs(E)<=100)+al.measured_offset_pixels(t),t));end
assert(all(beam>0),'qe_align_sequence_frames:Beam','Nonpositive beam reference');
elastic=squeeze(sum(F(abs(E)<=opts.q_centroid_energy,:,:),1));
central=abs(Q)<=opts.q_centroid_range+1e-12;
mu=sum(elastic(:,central).*Q(central),2)./sum(elastic(:,central),2);
drift=movmedian(mu-median(mu),opts.q_median_frames).';
if ~opts.do_q_drift,drift(:)=0;end
profile=sum(elastic,1);[~,p]=max(profile);jj=p+(-2:2);
w=profile(jj)/sum(profile(jj));x=Q(jj);probe=sqrt(max(sum(w.*(x-sum(w.*x)).^2)-dq^2/12,0));
A=struct('source',source,'raw_size',raw_size,'E',E,'Q',Q,'selected',selected,'F',F, ...
    'alignment',al,'beam',beam,'centroid_Ainv',mu,'drift_Ainv',drift,'sigma_probe',probe, ...
    'energy_origin_meV',Eraw(al.reference_pixel));
end
