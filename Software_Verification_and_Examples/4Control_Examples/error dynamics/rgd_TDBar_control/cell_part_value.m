function A = cell_part_value(A,block1,row_or_col,target_index)
a = block1(1);
b = block1(2);
A_temp = A{a,b};

if target_index == 2
    A_temp(:,row_or_col) = [];
else 
    A_temp(row_or_col,:) = [];
end
A{a,b} = A_temp;
end
