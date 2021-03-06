motdict = cell(0);
kek = unique(dump(:,1));
for i=1:length(kek)
    motdict{i,1} = kek{i};
    hurr = dump(ismember(dump(:,1),kek{i}),:);
    motdict{i,2} = '';
    wololo = '';
    motdict{i,3} = '';
    for j=1:size(hurr,1)
        holder = textscan(hurr{j,2},'%s','delimiter','-');
        motdict{i,2} = strcat(motdict{i,2},wololo,holder{1}{1});
        motdict{i,3} = strcat(motdict{i,3},wololo,hurr{j,3});
        wololo = ',';
    end
end
fid = fopen('motdict-parsed.txt','w');
for i=1:size(motdict,1)
    fprintf(fid,motdict{i,1});
    for j=2:size(motdict,2)
        fprintf(fid,['\t',motdict{i,j}]);
    end
    if i<size(motdict,1)
        fprintf(fid,'\n');
    end
end
fclose(fid);