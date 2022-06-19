function tenseg_tecplot(C,n_t,t_t,radius)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% input:connection matrix C, nodal coordinate n_t, members' force t_t,radius
% output:"tecplot_data.dat" file
% N node number,E element number,n node coord 3xN,e elem direction 3xE
%%
[E,N]=size(C);
T=size(n_t,2);
for i=1:E
    CC(i,:)=find(C(i,:));
end

fid=fopen('tecplot_data.dat','w');
for i=1:T
    n=reshape(n_t(:,i),[3,N]);
    e=n(:,CC(:,1))-n(:,CC(:,2));
    for j=1:E
        [~,~,v]=svd(e(:,j)');
        for I=0:1
        for J=0:1
        for K=0:1
            nn(:,8*(j-1)+4*I+2*J+K+1)=n(:,CC(j,I+1))+(-1)^J*v(:,2+K)*radius(j);
        end
        end
        end
    end
    NN=size(nn,2);
    fprintf(fid,'Title="tecplot data for time step %d"\n',i);
    fprintf(fid,'variables ="x","y","z","f","d"\n');
    fprintf(fid,'ZONE N =%d,E =%d,datapacking=block\n',NN,E);
    fprintf(fid,'varlocation=([4,6]=cellcentered)\n');
    fprintf(fid,'zonetype=febrick\n');
    
    fprintf(fid,'%7.5f\t',nn(1,:));fprintf(fid,'\n');
    fprintf(fid,'%7.5f\t',nn(2,:));fprintf(fid,'\n');
    fprintf(fid,'%7.5f\t',nn(3,:));fprintf(fid,'\n');
    fprintf(fid,'%15.5f\t',t_t(:,i));fprintf(fid,'\n');
    
    d=reshape(n_t(:,i)-n_t(:,1),[3,N]);
    dd=kron(sum(d(:,CC').^2).^0.5,ones(1,4));
    fprintf(fid,'%7.5f\t',dd);fprintf(fid,'\n');
    fprintf(fid,'%d %d %d %d\n',1:NN);
end
fclose(fid); 

end