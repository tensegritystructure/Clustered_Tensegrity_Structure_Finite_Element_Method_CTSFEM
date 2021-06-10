function [N_new,Cb_new,Cs_new] = catch_error(N,Cb,Cs)
% if iscell(Yt):
%     try Yt{1,1}*Cb
%     catch ME
%         if isa(ME,'double')
%             Yt
%         if sum(ismember('dimensions',ME.message)) == length(ismember('dimensions',ME.message))
%             try
try N*Cb';
catch ME
    if sum(ismember('dimensions',ME.message)) == length(ismember('dimensions',ME.message))
        N = N';
    end
end
N_new = N;
Cb_new = Cb;
Cs_new = Cs;
end