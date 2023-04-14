
function [energy, maxProbFold] = ProteinVQEObjective(parameters,hyperParams)    
% Construct and simulate the variational circuit 
ansatz = ProteinConfigAnsatz(parameters);
qState = simulate(ansatz);
allProbs = (qState.Amplitudes).^2;

% There are 10 qubits in the circuit, but only these 9 are used to define a fold 
foldQubits = [1:8 10];

% Get the most probable fold
% Only compute this if second output maxProbFold is requested
if nargout > 1 
    [~,idx] = max(allProbs);
    maxProbKet = char(qState.BasisStates(idx));
    maxProbFold = maxProbKet(foldQubits);
end

% Sample, and query/get the states and probabilities of the fold qubits 
qMeasurement = randsample(qState, hyperParams.numShots);
[states, probs] = querystates(qMeasurement, foldQubits);

% Sort the probabilities by the energy
[energies,sort_idx] = sort(exactHamiltonian(char(states), hyperParams));
probs = probs(sort_idx);

% Compute CVaR over the low energy tail of the energy distribution,
% delimited by a cutoff parameter alpha. 
alpha = .025;  
cut_idx = nnz(cumsum(probs) < alpha);
cvar_probs = probs(1:cut_idx);
cvar_probs(end+1) = alpha - sum(cvar_probs);

% Compute expectation energy as the sum of cutoff state energies weighted by their probability 
energy = dot(cvar_probs, energies(1:cut_idx+1))/alpha;
end