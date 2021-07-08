function Gp=reduced_control_variable(gr2,gr1)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the transforming matrix from two cell
%
% Inputs:
%   gr2: number corresponding to column
%   gr1: number corresponding to row
% Outputs:
%	Gp: a transforming matrix
gr1;gr2;
num=0;
for i=1:numel(gr2)
    num=num+numel(gr2{i});
end

Gp=eye(numel)
Gp=eye(size(C,1));
E=eye(size(C,1));
Gp1=[];
%% method 1

if ~isempty(gr)
   for i=1:numel(gr)       % this is to combine members in one group
Gp(:,gr{i}(1))=sum(E(:,gr{i}),2);
   end
      for i=1:numel(gr)         % set duplicate columns into 0
Gp(:,gr{i}(2:end))=0*Gp(:,gr{i}(2:end));
   end
end
% delete zero column
Gp(:,~max(Gp))=[];

end

