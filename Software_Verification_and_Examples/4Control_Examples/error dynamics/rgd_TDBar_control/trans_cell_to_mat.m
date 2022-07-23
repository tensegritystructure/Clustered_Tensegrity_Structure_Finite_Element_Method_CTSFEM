function mat_A = trans_cell_to_mat(A)
[a,b] = size(A{1,1});
mat_A = zeros(a,b,length(A(1,:)));
for i =1:length(A(1,:))
    mat_A(:,:,i) = A{1,i};
end
    
    