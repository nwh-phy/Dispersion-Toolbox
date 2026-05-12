function evidence = b1_peak_evidence_audit_classify(evidence, options)
%B1_PEAK_EVIDENCE_AUDIT_CLASSIFY Classify support for selected B1 peaks.
%
% The classifier is intentionally diagnostic: it labels whether a selected
% Lorentz branch point is directly supported by the local spectrum, only
% helped by tracking continuity, or suspicious enough to revisit.

arguments
    evidence table
    options.DataSupportScore (1,1) double = 0.70
    options.TrackingSupportScore (1,1) double = 0.30
    options.StrongComponentSSE (1,1) double = 0.10
    options.ModestComponentSSE (1,1) double = 0.03
    options.LowerReleaseMaxMeV (1,1) double = 120
    options.UpperSmallQReleaseMaxMeV (1,1) double = 150
    options.UpperPlateauReleaseMaxMeV (1,1) double = 120
    options.SmallQAbsAinv (1,1) double = 0.06
    options.RobustnessMaxMeV (1,1) double = 150
    options.UpperMaxGammaOverE (1,1) double = 1.4
    options.UpperMaxGammaMeV (1,1) double = 1600
end

n = height(evidence);
evidence = local_ensure_numeric(evidence, 'q_abs_Ainv', NaN(n, 1));
evidence = local_ensure_text(evidence, 'branch_label', repmat({''}, n, 1));
evidence = local_ensure_numeric(evidence, 'local_support_score', NaN(n, 1));
evidence = local_ensure_numeric(evidence, ...
    'component_sse_increase_fraction', NaN(n, 1));
evidence = local_ensure_numeric(evidence, ...
    'constraint_release_delta_meV', NaN(n, 1));
evidence = local_ensure_numeric(evidence, ...
    'robustness_max_delta_meV', NaN(n, 1));
evidence = local_ensure_numeric(evidence, 'gamma_over_E', NaN(n, 1));
evidence = local_ensure_numeric(evidence, 'gamma_meV', NaN(n, 1));

classes = repmat({'suspicious'}, n, 1);
reasons = repmat({''}, n, 1);
is_data = false(n, 1);
is_tracking = false(n, 1);
is_suspicious = true(n, 1);

for i = 1:n
    branch = lower(char(string(evidence.branch_label{i})));
    is_upper = contains(branch, 'upper');
    release_limit = options.LowerReleaseMaxMeV;
    if is_upper && evidence.q_abs_Ainv(i) <= options.SmallQAbsAinv
        release_limit = options.UpperSmallQReleaseMaxMeV;
    elseif is_upper
        release_limit = options.UpperPlateauReleaseMaxMeV;
    end

    release_delta = evidence.constraint_release_delta_meV(i);
    robust_delta = evidence.robustness_max_delta_meV(i);
    release_stable = isfinite(release_delta) && release_delta <= release_limit;
    robust_stable = ~isfinite(robust_delta) || ...
        robust_delta <= options.RobustnessMaxMeV;
    strong_component = evidence.component_sse_increase_fraction(i) >= ...
        options.StrongComponentSSE;
    modest_component = evidence.component_sse_increase_fraction(i) >= ...
        options.ModestComponentSSE;
    direct_support = evidence.local_support_score(i) >= ...
        options.DataSupportScore;
    weak_support = evidence.local_support_score(i) >= ...
        options.TrackingSupportScore;
    overbroad = is_upper && ((isfinite(evidence.gamma_over_E(i)) && ...
        evidence.gamma_over_E(i) > options.UpperMaxGammaOverE) || ...
        (isfinite(evidence.gamma_meV(i)) && ...
        evidence.gamma_meV(i) > options.UpperMaxGammaMeV));

    reason_parts = {};
    if ~release_stable
        reason_parts{end + 1} = 'release_unstable'; %#ok<AGROW>
    end
    if ~robust_stable
        reason_parts{end + 1} = 'robustness_unstable'; %#ok<AGROW>
    end
    if overbroad
        reason_parts{end + 1} = 'overbroad_upper_peak'; %#ok<AGROW>
    end
    if ~weak_support
        reason_parts{end + 1} = 'weak_local_support'; %#ok<AGROW>
    end
    if ~modest_component
        reason_parts{end + 1} = 'component_not_needed'; %#ok<AGROW>
    end

    if direct_support && strong_component && release_stable && ...
            robust_stable && ~overbroad
        classes{i} = 'data_supported';
        reasons{i} = 'local_support,strong_component,release_stable';
        is_data(i) = true;
        is_suspicious(i) = false;
    elseif weak_support && modest_component && release_stable && ...
            robust_stable && ~overbroad
        classes{i} = 'tracking_assisted';
        reasons{i} = 'weak_local_support,modest_component,release_stable';
        is_tracking(i) = true;
        is_suspicious(i) = false;
    else
        classes{i} = 'suspicious';
        if isempty(reason_parts)
            reason_parts = {'insufficient_combined_evidence'};
        end
        reasons{i} = strjoin(reason_parts, ',');
        is_suspicious(i) = true;
    end
end

evidence.evidence_class = classes;
evidence.evidence_reason = reasons;
evidence.is_data_supported = is_data;
evidence.is_tracking_assisted = is_tracking;
evidence.is_suspicious = is_suspicious;
end


function tbl = local_ensure_numeric(tbl, name, default_value)
if ~ismember(name, tbl.Properties.VariableNames)
    tbl.(name) = default_value;
    return
end
value = tbl.(name);
if iscell(value) || isstring(value) || ischar(value)
    tbl.(name) = str2double(string(value));
else
    tbl.(name) = double(value);
end
end


function tbl = local_ensure_text(tbl, name, default_value)
if ~ismember(name, tbl.Properties.VariableNames)
    tbl.(name) = default_value;
    return
end
value = tbl.(name);
if isstring(value)
    tbl.(name) = cellstr(value);
elseif ischar(value)
    tbl.(name) = cellstr(value);
end
end
