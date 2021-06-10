function [fig_out] = tenseg_plot_target_3D( N,C_b,C_s,target,fig_handle,highlight_nodes,view_vec )
% [fig_out] = TENSEG_PLOT( N,C_b,C_s,fig_handle,highlight_nodes,view_vec )
% creates a rough visualization figure for a given tensegrity structure
%
% Inputs:
%	N: node matrix (3 x n array for n nodes)
%	C_b (optional): bar connectivity matrix (beta x n array for beta bars)
%	C_s: string connectivity matrix (alpha x n array for alpha strings)
%	fig_handle (optional): figure in which to plot
%	highlight_nodes (optional): node(s) to highlight when plotting (vector
%		containing node numbers)
%	view_vec (optional): plotting view direction (see view())
%
% Outputs:
%	fig_out: figure handle in which plot was made
%
% Example:
%	tenseg_plot(N,C_b,C_b)

% Handle optional arguments
switch nargin
    case 6
        view_vec = [];
    case 5
        view_vec = [];
        highlight_nodes = [];
    case 4
        view_vec = [];
        highlight_nodes = [];
        fig_handle = [];
end
% Open specified figure or create new one
if isempty(fig_handle),
    fig_out = figure;
else
    fig_out = figure(fig_handle);
end
% y = zeros(length(target(:,1)),1);
% target_new = [target(:,1),y,target(:,2)];
target_new = target;

% Plot bar member vectors
if ~isempty(C_b),
    B = N*C_b';
    bar_start_nodes = zeros(3,size(B,2));
    for j = 1:size(B,2),
        bar_start_nodes(:,j) = N(:,C_b(j,:)==-1);
    end
    % 	quiver3(bar_start_nodes(1,:),bar_start_nodes(2,:),bar_start_nodes(3,:),B(1,:),B(2,:),B(3,:),'black.','Autoscale','off','LineWidth',2)
    quiver3(bar_start_nodes(1,:),bar_start_nodes(2,:),bar_start_nodes(3,:),B(1,:),B(2,:),B(3,:),'color','[0.5 0.5 0.5].','Autoscale','off','LineWidth',3,'ShowArrowHead','off');
    hold on
end
% Add axes labels
xlabel('x'); ylabel('y'); zlabel('z');
% Plot string member vectors
if ~isempty(C_s),
    S = N*C_s';
    string_start_nodes = zeros(3,size(S,2));
    for j = 1:size(S,2),
        string_start_nodes(:,j) = N(:,C_s(j,:)==-1);
    end
    % 	quiver3(string_start_nodes(1,:),string_start_nodes(2,:),string_start_nodes(3,:),S(1,:),S(2,:),S(3,:),'red.','Autoscale','off');
    quiver3(string_start_nodes(1,:),string_start_nodes(2,:),string_start_nodes(3,:),S(1,:),S(2,:),S(3,:),'color','[1 0.7529 0.7961].','Autoscale','off','LineWidth',2,'ShowArrowHead','off');
    hold on
end

% Plot node markers
% plot3(N(1,:),N(2,:),N(3,:),'black.','MarkerSize',6)
% plot3(N(1,:),N(2,:),N(3,:),'color',[1 0.7529 0.7961],'MarkerSize',6)

% Highlight specified nodes if applicable
for j=1:numel(highlight_nodes),
    node_index = highlight_nodes(j);
    plot3(N(1,node_index),N(2,node_index),N(3,node_index),'kd','MarkerSize',4,'MarkerFaceColor','blue')
end
plot3(target_new(:,1),target_new(:,2),target_new(:,3),'bo','MarkerSize',4,'MarkerFaceColor','blue')
% Modify plot display
grid off; axis equal;
if isempty(view_vec),
    [~, view_vec_derived] = tenseg_axisview(N);
    view_vec = view_vec_derived;
end
view(view_vec)

end