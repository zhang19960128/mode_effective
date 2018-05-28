borne=importdata("born.txt");
freq=importdata("fre.txt");
len=length(freq);
epsall=zeros(3,3);
fileID = fopen('alldie.txt','w');
for i=4:len
    add=die(borne(i,1:end),freq(i));
    add
    for j=1:9
        fprintf(fileID,'%f ',add(j));
    end
    fprintf(fileID,'\n');
    epsall=epsall+add;
end
epsall
fclose(fileID);