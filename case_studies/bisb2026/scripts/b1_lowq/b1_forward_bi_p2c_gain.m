function [X,lim] = b1_forward_bi_p2c_gain(B,mode)
%B1_FORWARD_BI_P2C_GAIN  Smooth odd exit response, or a derivative shift proxy.
% X times its coefficient gives log(measured even-signal multiplier).
E = B.E;
q = B.q(:).';
u = abs(q)/max(abs(q));
switch mode
    case 'none'
        X = zeros(numel(E),numel(q),size(B.Y,3),0);
        lim = [];
    case 'gainE'
        z = log(max(E,250)/700);
        k = sign(q).*[u; u.*(1-u); u.*(1-u).^2]/2;
        k = reshape(k.',1,numel(q),1,3);
        X = cat(4,k.*ones(size(E)),k.*z);
        X = repmat(X,1,1,size(B.Y,3),1);
        lim = [1.2 4 8 1.2 4 8];
    case 'shift'
        jp = find(q>0);
        [tf,jm] = ismembertol(-q(jp),q,1e-10,'DataScale',1);
        assert(all(tf));
        m = movmean(log(max((B.archive(:,jp,:)+B.archive(:,jm,:))/2,1)),9,1);
        d = zeros(size(m));
        % Hold the derivative shape at the B1 lower edge for extrapolation.
        for i=1:numel(jp)
            a = max(i-1,1); z = min(i+1,numel(jp));
            d(:,i,:) = (m(:,z,:)-m(:,a,:))/(q(jp(z))-q(jp(a)));
        end
        k = find(E>=250,1);
        d(1:k-1,:,:) = repmat(d(k,:,:),k-1,1,1);
        v = zeros(size(B.Y)); v(:,jp,:) = -d; v(:,jm,:) = d;
        X = cat(4,v*1e-3,v.*(max(E,250)*B.intr.qE_per_meV));
        lim = [.6 1];
    otherwise
        error('Unknown odd response %s.',mode);
end
end
