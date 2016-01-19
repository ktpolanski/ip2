fid = fopen('/Users/krzysztofpolanski/Downloads/Arabidopsis_thaliana_2015_05_29_7-38_am/TF_Information_all_motifs.txt','r');
fgetl(fid);
dump = cell(0);
while ~feof(fid)
    line = fgetl(fid);
    line = textscan(line,'%s','delimiter','\t');
    dump(end+1,1:length(line{1}))=line{1};
end
fclose(fid);
dump = dump(:,[4 6 7]);

fid = fopen('/Users/krzysztofpolanski/Documents/python/ip2/HMT/murray_genes/individual_gene_motif_hits/all_motifs.txt','r');
mots = cell(0);
while ~feof(fid)
    mots{end+1} = fgetl(fid);
end
fclose(fid);

mask = ismember(dump(:,1),mots);
dump = dump(mask,:);
[~,inds] = sort(dump(:,1));
dump = dump(inds,:);

fid = fopen('/Users/krzysztofpolanski/Documents/python/ip2/HMT/murray_genes/motif_dictionary.txt','w');
for i=1:size(dump,1)
    fprintf(fid,dump{i,1});
    for j=2:size(dump,2)
        fprintf(fid,['\t',dump{i,j}]);
    end
    if i<size(dump,1)
        fprintf(fid,'\n');
    end
end
fclose(fid);