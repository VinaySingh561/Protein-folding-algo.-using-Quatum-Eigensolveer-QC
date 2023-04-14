function mat = buildMJInteractions(protein)
% Create the MJ interaction energy matrix for the protein, specified with
% 1-letter codes 
N = length(protein);
mat = zeros(N,N);
rng(29507)
MJ = rand(20)*-6;
MJ = triu(MJ) + triu(MJ, 1)';
acids = ["C","M","F","I","L","V","W","Y","A","G","T","S","N","Q","D","E","H","R","K","P"];
acid2idx = dictionary(acids, 1:20);
for i = 1:N
    for j=1:N
        mat(i,j) = MJ(acid2idx(protein(i)), acid2idx(protein(j)));
    end
end
end

