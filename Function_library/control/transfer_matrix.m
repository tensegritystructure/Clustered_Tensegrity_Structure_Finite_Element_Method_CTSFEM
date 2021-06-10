function Ic=transfer_matrix(b,a)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% Inputs:
% a:vector with indexs (longer)
% b:vector with indexs (shorter)
% Output:
% Ic: is the transfer matrix 
% Satisfy    b=Ic'*a

Ic=zeros(numel(a),numel(b));
for i=1:numel(b)
    term=find(a==b(i));
    Ic(term,i)=1;
end

