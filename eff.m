data_vector=importdata("modeone.txt");
data_born=importdata("born.txt");
atom_num=length(data_born)/3;
atom_num=int64(atom_num);
atom_mass=[];
for i=1:int64(atom_num/5*3)
    atom_mass=[atom_mass,15.999];
end
for i=1:int64(atom_num/5*1)
    atom_mass=[atom_mass,91.224];
end
for i=1:int64(atom_num/5*1)
    atom_mass=[atom_mass,137.33];
end
z_star=zeros(1,3);
addall=[];
for k=1:3
for i=1:atom_num
    i_born_tensor=data_born(int64(3*(i-1)+1):int64(3*(i)),1:end);
    i_vector=data_vector(i,1:end);
    for j=1:3
        add=i_born_tensor(k,j)*i_vector(j)/sqrt(atom_mass(i));
        z_star(k)=z_star(k)+add;
        if(k==1)
            addall=[addall;add]
        end
    end
end
end
z_star
