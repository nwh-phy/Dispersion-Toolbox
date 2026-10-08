function bins = qe_prepare_count_bins(qe, N, options)
%QE_PREPARE_COUNT_BINS Fixed nonoverlapping native-channel bins.
% Missing samples exclude a whole q channel (common energy support). No
% omitnan averaging, zero padding, cross-side or across-gap combinations.
% Variance propagation is conditional on independent input channels; leave
% it NaN when unknown. Between-channel scatter is NOT measurement variance.
arguments
    qe struct
    N (1,1) double {mustBePositive,mustBeInteger}
    options.q_range_Ainv (1,2) double = [-0.015 0.015]
    options.q_skip_Ainv (1,1) double {mustBeNonnegative} = 0.0005
    options.invalid_q logical = false(1,size(qe.intensity,2))
    options.variance double = []
    options.variance_source string = "unknown"
    options.processing_level string = "unresolved"
end
q=double(qe.q_Ainv(:).'); X=double(qe.intensity); dq=qe.dq_Ainv;
assert(size(X,2)==numel(q) && all(isfinite(q)) && all(diff(q)>0), ...
    'qe_prepare_count_bins:Axis','q must be finite, ascending and match columns');
assert(isfinite(dq)&&dq>0,'qe_prepare_count_bins:Dq','Invalid dq');
assert(numel(options.invalid_q)==numel(q),'qe_prepare_count_bins:Mask','Mask size mismatch');
source=1:numel(q);
if isfield(qe,'source_channel'), source=qe.source_channel(:).'; end
assert(numel(source)==numel(q)&&all(diff(source)>0),'qe_prepare_count_bins:Channels','Invalid native channels');
V=options.variance;
if isempty(V), V=nan(size(X)); else
 assert(isequal(size(V),size(X)) && all(V(isfinite(V))>=0),'qe_prepare_count_bins:Variance','Invalid variance');
end
tol=dq*1e-7;
valid=q>=options.q_range_Ainv(1)-tol & q<=options.q_range_Ainv(2)+tol & ...
 abs(q)>=options.q_skip_Ainv-tol & q~=0 & ~options.invalid_q(:).' & all(isfinite(X),1);
idx=find(valid); groups={}; run=[];
for j=idx
 if ~isempty(run) && (j~=run(end)+1 || source(j)~=source(run(end))+1 || ...
        sign(q(j))~=sign(q(run(end))) || abs(q(j)-q(run(end))-dq)>tol)
  groups=append(groups,run,N); run=[];
 end
 run(end+1)=j; %#ok<AGROW>
end
groups=append(groups,run,N);
bins=struct('units',struct([]),'sum',zeros(size(X,1),0),'mean',zeros(size(X,1),0), ...
 'variance_sum',zeros(size(X,1),0),'variance_mean',zeros(size(X,1),0), ...
 'member_scatter_variance',zeros(size(X,1),0),'valid_count',zeros(size(X,1),0), ...
 'valid_q',valid,'variance_source',options.variance_source, ...
 'processing_level',options.processing_level,'energy_meV',qe.energy_meV);
for k=1:numel(groups)
 ii=groups{k}; n=numel(ii); center=mean(q(ii));
 u=struct('bin_id',k,'q_Ainv',center,'q_abs_Ainv',abs(center),'q_indices',ii, ...
  'source_mode',sprintf('combined_q_binning_%d',N),'source_q_count',n, ...
  'source_q_Ainv',char(strjoin(compose('%.12g',q(ii)),',')), ...
  'source_q_index',char(strjoin(compose('%d',ii),',')), ...
  'source_channel',source(ii),'bin_size_requested',N,'partial',n<N, ...
  'q_left',q(ii(1))-dq/2,'q_right',q(ii(end))+dq/2,'q_width',n*dq);
 if N==1, u.source_mode='single_q_direct'; end
 if k==1, bins.units=u; else, bins.units(k)=u; end
 bins.sum(:,k)=sum(X(:,ii),2); bins.mean(:,k)=bins.sum(:,k)/n;
 bins.variance_sum(:,k)=sum(V(:,ii),2); bins.variance_mean(:,k)=bins.variance_sum(:,k)/n^2;
 bins.member_scatter_variance(:,k)=var(X(:,ii),0,2);
 bins.valid_count(:,k)=n;
end
end
function groups=append(groups,run,N)
for j=1:N:numel(run), groups{end+1}=run(j:min(j+N-1,end)); end
end
