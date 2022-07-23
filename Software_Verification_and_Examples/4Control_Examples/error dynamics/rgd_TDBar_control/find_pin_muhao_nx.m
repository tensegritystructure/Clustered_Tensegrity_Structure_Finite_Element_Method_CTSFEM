function pinned_nodes = find_pin_muhao_nx(N,x)
pinned_nodes = [];
for i = 1:length(N(1,:))
    if (N(1,i) <x)
%         &(N(3,i)==0)
        pinned_nodes = [pinned_nodes i];
    end
end