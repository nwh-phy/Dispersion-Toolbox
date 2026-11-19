function [c,o,hScan] = b1_forward_bi_p2b_config(id)
%B1_FORWARD_BI_P2B_CONFIG  Registered geometry, archived P2a response and B1 starts.
a = b1_forward_bi_p2a_config(id);
p = a.response.paths;
tag = b1_sha256(a.response);
p.response_dir = fullfile(p.input_dir,['response_' tag(1:8)]);
d = bisb_datasets(id);
nodes = struct('E0',[d.dq .001 .0015 .002 .0025 .0035 .005], ...
    'G',[d.dq .002 .005], 'A',[d.dq .002 .005]);
c = b1_p2b_config(p,nodes);
o = b1_p2a_operator_config(p);
c.highq.file = '';
c.lobo_block = 6;
c.predict_block = true;
c.version = 'P2b-Bi-2026-11-19';
c.Q_curve = unique([c.Q_curve c.nodes.E0]);
if o.h > .005
    hScan = [0.0033 0.00495 0.0066 0.00825 0.0099];
else
    hScan = [0.00165 0.002475 0.0033 0.004125 0.00495 0.0066];
end
switch id
    case 'Bi1107_10w'
        e = [650 650 640 625 620 620 620;1400 1410 1425 1440 1450 1490 1525];
        g = [300 450 550;650 1200 1550];
    case 'Bi1107_20w'
        e = [620 620 620 620 620 620 620;1400 1420 1440 1460 1480 1510 1560];
        g = [420 550 650;850 1400 1750];
    case 'Bi1117_10w'
        e = [605 605 610 610 615 620 630;1350 1380 1400 1420 1440 1480 1520];
        g = [300 420 520;650 1150 1500];
    case 'Bi1117_20w'
        e = [585 590 595 600 605 610 620;1415 1430 1445 1460 1475 1500 1540];
        g = [300 420 520;650 1150 1500];
    otherwise
        error('b1_forward_bi_p2b_config:Session','Unsupported acquisition %s.',id);
end
flat = repmat(median(e,2),1,size(e,2));
if strcmp(id,'Bi1107_20w'), flat = repmat([700;1650],1,size(e,2)); end
drop = e;
drop(1,1:3) = [450 550 (550+e(1,4))/2];
c.starts = struct('name',{'P2a handback','flat','low-Q drop'},'E0',{e,flat,drop},'G',{g,g,g});
end
