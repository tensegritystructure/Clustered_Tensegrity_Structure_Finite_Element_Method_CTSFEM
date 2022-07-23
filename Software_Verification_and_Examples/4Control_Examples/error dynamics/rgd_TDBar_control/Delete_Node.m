function [N,CBT,CST]=Delete_Node(delete_index,N,CBT,CST)
% [N,CBT,CST]=Delete_Node(delete_index,N,CBT,CST)  helps user to delete a 
% specific node and all the bars and strings connected to that node.

% Inputs:
%	delete_index: specify which node to be deleted
%   N: Nodal Matrix 
%   CBT: Bar Connectivity Matrix
%   CST: String Connectivity Matrix
% Outputs:
%   N: New Nodal Matrix 
%   CBT: New Bar Connectivity Matrix
%   CST: New String Connectivity Matrix

%% Dimension correction
CBT = CBT';CST = CST';
%% Remove the bars and strings connect to the specific node.
CBT(:,find(abs(CBT(delete_index,:))==1)) = [];
CST(:,find(abs(CST(delete_index,:))==1)) = [];
%% Remove the node appereas in row of CBT anc CST 
N(:,delete_index)= []; CBT(delete_index,:) = []; CST(delete_index,:) = [];
%% Dimension correction
CBT = CBT';CST = CST';

end