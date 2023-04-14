function energies = exactHamiltonian(bitstrings,hyperParams)
% Compute the Hamiltonian for each bit string (i.e., the energy for each fold).
% See [2] for details. This does not consider the Hch constraint from
% side-chains and the interaction term is only 1-nearest-neighbor (1-NN).

% H = Hgc + Hin_1

lambdaDis = 720;    % Penalty for interaction distance
lambdaLoc = 20;     % Penalty for interaction location
lambdaBack = 50;    % Penalty for unphysical geometry

energies = zeros(size(bitstrings,1),1);
numBeads = length(hyperParams.protein);

for idx = 1:length(energies)
    bitstring = bitstrings(idx,:);
    config = hyperParams.turn2qubit; 
    config(config=='q') = bitstring(1:hyperParams.numQubitsConfig);
    turns = bin2dec(reshape(config,2,numBeads-1)');
    
    %% Geometric Hamiltonian Hgc

    % Count number of adjacent turns which are equal and impose a penalty for each
    energies(idx) = lambdaBack*sum(turns(1:end-1) == turns(2:end));
  
    %% 1-NN Interaction Hamiltonian Hin
    currInteractionQubit = hyperParams.numQubitsConfig;
    for i=1:(numBeads-4) 
       for j=(i+5):2:numBeads 
           
           currInteractionQubit = currInteractionQubit+1;
           if bitstring(currInteractionQubit)=='0'
               continue;
           end       
            
           % Add the interaction energy 
           energies(idx) = energies(idx) + hyperParams.interactionEnergy(i,j);
            
           % Compute distances between interacting beads
           deltaN_ij = zeros(1,4);
           deltaN_ir = zeros(1,4);
           deltaN_mj = zeros(1,4);
           for k=0:3
               deltaN_ij(k+1) = (-1).^(i:j-1)*(turns(i:j-1) == k);
               deltaN_ir(k+1) = (-1).^(i:j-2)*(turns(i:j-2) == k);
               deltaN_mj(k+1) = (-1).^(i+1:j-1)*(turns(i+1:j-1) == k);
           end
           d_ij = norm(deltaN_ij)^2;
           d_ir = norm(deltaN_ir)^2;
           d_mj = norm(deltaN_mj)^2;
            
           % Add penalty for distance not equal to 1
           energies(idx) = energies(idx) + lambdaDis*(d_ij-1);         
            
           % Add penalty for unphysical nearest neighbour collisons
           energies(idx) = energies(idx) + lambdaLoc*(2-d_ir);
           energies(idx) = energies(idx) + lambdaLoc*(2-d_mj);
           if i-1 >= 1
               for k=0:3
                   deltaN_mj(k+1) = (-1).^(i-1:j-1)*(turns(i-1:j-1) == k);
               end
               d_mj = norm(deltaN_mj)^2;
               energies(idx) = energies(idx) + lambdaLoc*(2-d_mj);
           end
           
           if j+1 <= numBeads
               for k=0:3
                   deltaN_ir(k+1) = (-1).^(i:j)*(turns(i:j) == k);
               end
               d_ir = norm(deltaN_ir)^2;
               energies(idx) = energies(idx) + lambdaLoc*(2-d_ir);
           end
       end
    end
end
end

