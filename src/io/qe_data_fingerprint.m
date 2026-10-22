function hash = qe_data_fingerprint(qe)
%QE_DATA_FINGERPRINT Fingerprint projected counts and their calibrated axes.
md = java.security.MessageDigest.getInstance('SHA-256');
fields = {'intensity', 'energy_meV', 'q_Ainv', 'q_channel', 'q_zero_index', 'dq_Ainv'};
for k = 1:numel(fields)
    values = double(qe.(fields{k}));
    md.update(typecast(uint64(size(values)), 'int8'));
    % Bound the temporary typecast allocation for large q-E matrices.
    for first = 1:524288:numel(values)
        last = min(first + 524287, numel(values));
        md.update(typecast(values(first:last), 'int8'));
    end
end
hash = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
