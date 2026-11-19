function [A, q] = b1_forward_stage_frames(file, C, qmax)
%B1_FORWARD_STAGE_FRAMES  Restore a native E x q x frame cube including q=0.
m=matfile(file);
ix=find(abs(C.q)<=qmax+1e-9);
q=[C.q(ix),0]; [q,order]=sort(q);
F=cat(3,m.frames(:,:,ix),m.center_frames);
A=permute(double(F(:,:,order)),[1 3 2]);
end
