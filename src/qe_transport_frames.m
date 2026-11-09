function [frames, max_weight_error] = qe_transport_frames(F, Q, dq, centers, widths, drift, expected)
%QE_TRANSPORT_FRAMES  Bin aligned frames in q while following the beam drift of each frame.
%   frames = qe_transport_frames(F, Q, dq, centers, widths, drift) returns
%   E x frame x bin counts (single). Bin b of frame t collects the native
%   channels Q (width dq) weighted by their overlap with
%       [centers(b) + drift(t) - widths(b)/2, centers(b) + drift(t) + widths(b)/2],
%   so each bin keeps the same physical q while the beam moves. F is
%   E x frame x channel (e.g. qe_align_sequence_frames). For native channels use
%   centers = Q and widths = dq.
%   expected (default widths/dq) is the number of channels a bin should collect;
%   max_weight_error = max |sum of weights - expected| (counts are conserved
%   when it is ~0; callers usually assert it below 1e-8).
%
%   See also: qe_align_sequence_frames, qe_prepare_count_bins
arguments
    F {mustBeNumeric}
    Q (1,:) double
    dq (1,1) double {mustBePositive}
    centers (1,:) double
    widths (1,:) double
    drift (1,:) double
    expected (1,:) double = widths / dq
end
nE=size(F,1);T=size(F,2);B=numel(centers);
assert(numel(widths)==B&&numel(expected)==B&&numel(drift)==T&&size(F,3)==numel(Q), ...
    'qe_transport_frames:Size','Inconsistent sizes of F, Q, bins and drift');
frames=zeros(nE,T,B,'single');max_weight_error=0;
for b=1:B
    for t=1:T
        left=centers(b)+drift(t)-widths(b)/2;right=centers(b)+drift(t)+widths(b)/2;
        w=max(0,min(Q+dq/2,right)-max(Q-dq/2,left))/dq;
        max_weight_error=max(max_weight_error,abs(sum(w)-expected(b)));
        used=find(w>1e-12);frames(:,t,b)=single(reshape(F(:,t,used),nE,[])*w(used).');
    end
end
end
