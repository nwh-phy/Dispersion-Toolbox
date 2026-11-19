function [models,profile,limits,base] = b1_forward_bi_p2c_search(B,P,held,seed,doLimits,depth)
%B1_FORWARD_BI_P2C_SEARCH  Matched material families and weight profiles.
% Grid model shapes are compared with independent nuisance refits. Weight
% and weight-ratio constraints are imposed on the linear design directly.
t = tic;
if nargin<6, depth=0; end
ne=numel(B.cfg.nodes.E0);
nj=B.intr.n_branch;
ng=numel(B.cfg.nodes.G);
if P.fullB1
    free=1:nj*(ne+ng);
else
    free=[1 2 nj*ne+1 nj*ne+2];
    if strcmp(P.id,'Bi1107_20w') && nj>1
        free=[free nj*ne+(ng+1:2*ng)];
    end
end
C=struct('free',free);
[base,D]=b1_forward_bi_p2c_fit(B,P,held,seed,[],C);
C.design=D;
fl = C; fl.lite=true; fl.maxIter=24;
fixed=fl; fixed.free=[];
fixed.cols0=local_cols(B,base);
models=cell(3,1); models{1}=base;
profile={}; bounds={};
fams={'interband','drude'};
for k=1:numel(fams)
    fam=fams{k};
    [grid,width]=local_grid(P,fam);
    allRows=cell(numel(grid)*numel(width),7);
    i=0;
    for j=1:numel(grid)
        for w=width
            i=i+1;
            shape=[grid(j) w];
            [L,frac,normL]=b1_forward_bi_p2c_material(B,fam,shape,P.window(1));
            G=b1_forward_bi_p2c_fit(B,P,held,base,L,fixed);
            allRows(i,:)={shape,L,frac,G,false,j,normL};
        end
    end
    % Refit all threshold/relaxation grid points. Refine the second and
    % third widths if their fixed-shape fits are close to the best width.
    for j=1:numel(grid)
        group=find([allRows{:,6}]==j);
        frozen=cellfun(@(f) f.objective,allRows(group,4));
        [~,ord]=sort(frozen);
        near=frozen(ord)-frozen(ord(1))<=max(8,4*base.chi2dof);
        candidates=group(ord(near));
        for id=candidates(:).'
            G=b1_forward_bi_p2c_fit(B,P,held,base,allRows{id,2},fl);
            allRows{id,4}=G; allRows{id,5}=true;
        end
    end
    scored=find([allRows{:,5}]);
    vals=cellfun(@(f) f.objective,allRows(scored,4));
    [~,j]=min(vals);
    best=scored(j);
    model=b1_forward_bi_p2c_fit(B,P,held,allRows{best,4},allRows{best,2},C);
    model.shape=allRows{best,1}; model.commonFrac=allRows{best,3};
    model.amplitude=model.W/allRows{best,7};
    allRows{best,4}=model;
    models{k+1}=model;
    scale=max(1,model.chi2dof);
    for j=1:size(allRows,1)
        G=allRows{j,4}; x=allRows{j,1};
        profile(end+1,:)={fam,held,x(1),x(2),allRows{j,5},G.chi2,G.objective, ...
            G.objective-model.objective,(G.objective-model.objective)/scale,G.W,G.relativeW,G.flag}; %#ok<AGROW>
    end
    if doLimits
        vals=cellfun(@(f) f.objective,allRows(scored,4));
        active=scored(vals <= model.objective+2.705543*scale);
        target=model.objective+2.705543*scale;
        for field={'W','ratio'}
            metric=field{1};
            cand=zeros(numel(active),1); closed=false(size(cand));
            for jj=1:numel(active)
                id=active(jj); G=allRows{id,4};
                fixed.cols0=local_cols(B,G);
                [cand(jj),closed(jj)] = local_upper(B,P,held,G,allRows{id,2},fixed,metric,target,[],false);
            end
            factors=[allRows{active,3}].';
            provisional=[cand,cand.*factors];
            check=unique([local_top(provisional(:,1));local_top(provisional(:,2))]);
            for jj=check(:).'
                id=active(jj); G=allRows{id,4};
                [cand(jj),closed(jj),delta] = local_upper(B,P,held,G,allRows{id,2},fl,metric,target,cand(jj),true);
                s=allRows{id,1};
                bounds(end+1,:)={fam,held,metric,s(1),s(2),true,cand(jj),cand(jj)*factors(jj),closed(jj),delta,scale}; %#ok<AGROW>
            end
            remaining=setdiff(1:numel(active),check);
            for j=1:numel(remaining)
                jj=remaining(j);
                s=allRows{active(jj),1};
                bounds(end+1,:)={fam,held,metric,s(1),s(2),false,cand(jj),cand(jj)*factors(jj),closed(jj),NaN,scale}; %#ok<AGROW>
            end
        end
    end
    fprintf('%s held %d [%d,%d] %s %s: grid best [%g %g], W %.5g, train dchi2 %.4g, held %.5g, time %.1fs\n', ...
        P.id,held,P.window,P.tail,fam,model.shape,model.W,base.chi2-model.chi2,model.held_chi2n,toc(t));
end
profile=cell2table(profile,'VariableNames',{'family','held','parameter','width','refit','chi2','objective', ...
    'delta_objective','delta_scaled','W','W_over_low','flag'});
limits=cell2table(reshape(bounds,[],11),'VariableNames',{'family','held','metric','parameter','width', ...
    'refit','upper','upper_common','closed','bracket_delta','scale'});
if isfield(P,'nesting_check') && P.nesting_check
    challenger=base;
    for k=2:3
        F=b1_forward_bi_p2c_fit(B,P,held,models{k},[],C);
        if F.objective<challenger.objective, challenger=F; end
    end
    tol=max(.25,1e-5*base.objective);
    if challenger.objective<base.objective-tol && depth<2
        fprintf('%s held %d: cross-start lowers nested null objective by %.6g. Reprofiling.\n', ...
            P.id,held,base.objective-challenger.objective);
        [alternate,pa,ua,base]=b1_forward_bi_p2c_search(B,P,held,challenger,doLimits,depth+1);
        models{1}=base;
        for k=1:2
            if alternate{k+1}.objective<=models{k+1}.objective
                models{k+1}=alternate{k+1};
                take=strcmp(pa.family,fams{k});
                profile=[profile(~strcmp(profile.family,fams{k}),:);pa(take,:)];
                if height(ua)>0   % folds without limits return an empty table with no family names
                    take=strcmp(ua.family,fams{k});
                    limits=[limits(~strcmp(limits.family,fams{k}),:);ua(take,:)];
                end
            end
        end
    elseif challenger.objective<base.objective
        base=challenger; models{1}=base;
    end
end
end

function [grid,width]=local_grid(P,fam)
if isfield(P,'briefGrid') && P.briefGrid
    if strcmp(fam,'interband')
        grid=[0 150 200 250 300 400 550 800];
        width=[25 85];
    else
        grid=[25 60 120 240 500];
        width=NaN;
    end
elseif strcmp(fam,'interband')
    grid=P.thresholdGrid;
    width=[20 55 110];
else
    grid=[20 40 70 105 160 240 360 540];
    width=NaN;
end
end

function A=local_cols(B,F)
th=b1_p2b_decode(F.u,B.cfg.nodes,B.intr.n_branch);
C=qe_lowq_joint_forward('eval',B.cache,th,B.intr);
A=reshape(permute(C,[1 3 4 2]),[],size(C,2));
end

function ix=local_top(x)
[~,ix]=maxk(x,min(2,numel(x)));
end

function [hi,closed,delta]=local_upper(B,P,held,seed,L,C,metric,target,hint,full)
% Bracket the positive crossing using the penalized objective, with the
% nuisance free list fixed by C. Store the outer bracket endpoint.
if strcmp(metric,'W')
    w0=seed.W;
    oneCount=1/max(max(abs(L),[],'all'),1);
    units=max([seed.Wlow,abs(w0),oneCount]);
else
    if seed.Wlow<=0, hi=Inf; closed=false; delta=NaN; return; end
    w0=seed.relativeW;
    units=max(1,abs(w0));
end
units=max(units,realmin);
minStep=units*1e-12;
lo=max(w0,0);
if isempty(hint)
    step=max(units*1e-4,w0*.02);
    hi=lo+step;
else
    hi=max(hint,lo+max(abs(lo)*.001,minStep));
end
closed=false; delta=NaN;
Fprev=seed;
if full
    C.maxIter=18;
    it=16; expand=18; tol=.002;
else
    it=22; expand=24; tol=.001;
end
for i=1:expand
    [d,G]=obj(hi,Fprev);
    if d>=0, closed=true; break; end
    if full, Fprev=G; end
    step=max(hi-lo,minStep);
    lo=hi;
    hi=hi+2*step;
end
if ~closed, hi=Inf; return; end
delta=d;
if full, Fprev=G; end
for i=1:it
    mid=(lo+hi)/2;
    [d,G]=obj(mid,Fprev);
    if full, Fprev=G; end
    if d>=0
        hi=mid; delta=d;
    else
        lo=mid;
    end
    if hi-lo<tol*max(abs(hi),minStep), break; end
end
    function [d,G]=obj(x,start)
        C.(metric)=x;
        G=b1_forward_bi_p2c_fit(B,P,held,start,L,C);
        d=G.objective-target;
    end
end
