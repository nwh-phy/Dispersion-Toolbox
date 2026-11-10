function b1_fullq_run_part(run_root, session, part)
% Split the 78 bins of a session into three workers with similar load:
% 1 = q<0 with |q|<=0.0295 plus q>0 with |q|>=0.046, 2 = the mirror set,
% 3 = 0.031 <= |q| <= 0.0445 on both sides.
S=load(fullfile(run_root,'stage1',sprintf('frames_%s.mat',session)),'D'); q=S.D.q_center; a=abs(q); t=1e-9;
switch part
 case 1, sel=(q<0&a<=0.0295+t)|(q>0&a>=0.046-t);
 case 2, sel=(q>0&a<=0.0295+t)|(q<0&a>=0.046-t);
 case 3, sel=a>=0.031-t&a<=0.0445+t;
end
b1_fullq_fit_worker(run_root,session,find(sel));
end
