function dq_Ainv = infer_qe_dq_Ainv(source_path)
%INFER_QE_DQ_AINV  Infer BiSb q-EELS momentum step from session path.
%   The 20260120 BiSb sessions use 0.0005 1/A per pixel for 10w and
%   0.00025 1/A per pixel for 20w. Return NaN when the path is not one of
%   the calibrated session families.

lower_path = lower(char(source_path));

if contains(lower_path, '20w')
    dq_Ainv = 0.00025;
elseif contains(lower_path, '10w')
    dq_Ainv = 0.0005;
else
    dq_Ainv = NaN;
end
end
