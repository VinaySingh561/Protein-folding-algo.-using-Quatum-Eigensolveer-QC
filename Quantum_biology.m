%Configuration Qubits
hyperParams.protein = 'VINAYSINGH'; 
turn2qubit_prefix = '0100q1';

% Define the number of qubits to add to the turn2qubit string
num_qubits_to_add = length(hyperParams.protein) - 4;

% Define the string to be added to the turn2qubit string
qubit_string_to_add = repmat('qq', 1, num_qubits_to_add);

% Concatenate the prefix and the qubit string to create the final turn2qubit string
hyperParams.turn2qubit = [turn2qubit_prefix, qubit_string_to_add];
            
hyperParams.numQubitsConfig = sum(hyperParams.turn2qubit=='q');

% Interaction Qubits
hyperParams.numQubitsInteraction = 2;

% Call the buildMJInteractions function to find the interaction energies for the protein.
hyperParams.interactionEnergy = buildMJInteractions(hyperParams.protein);

% Compute Minimum Energy for All Folds
hyperParams.numQubitsTotal = hyperParams.numQubitsConfig + hyperParams.numQubitsInteraction; 

allFolds = dec2bin(0:2^hyperParams.numQubitsTotal-1,hyperParams.numQubitsTotal);
allEnergies = exactHamiltonian(allFolds,hyperParams);

hyperParams.GroundState.Energy = min(allEnergies);
hyperParams.GroundState.Index = find(allEnergies == hyperParams.GroundState.Energy);
allFolds(hyperParams.GroundState.Index,:)

% Find the energy of each of the identified folds.
allEnergies(hyperParams.GroundState.Index)

%Define the number of shots to use in ProteinVQEObjective, and create a function handle that passes in all of the parameter values to ProteinVQEObjective.
hyperParams.numShots = 1024; 
objFcn = @(theta) ProteinVQEObjective(theta,hyperParams);


% Call ProteinConfigAnsatz with random angles to construct the circuit. Plot the circuit to view the qubits and gates.
ansatz = ProteinConfigAnsatz(rand(2,9));
plot(ansatz)


numAngles = 2*hyperParams.numQubitsTotal; 
rng default

options = optimoptions("surrogateopt",...
    "MaxFunctionEvaluations",10, ...
    "PlotFcn","optimplotfval",...
    "InitialPoints",pi*ones(numAngles,1));

lb = repmat(-pi,numAngles,1);
ub = repmat(pi,numAngles,1);
[angles,minEnergy] = surrogateopt(objFcn,lb,ub,[],[],[],[],[],options);

[groundStateEnergy,groundStateFold] = ProteinVQEObjective(angles,hyperParams)

allFolds(allEnergies==minEnergy,:)

plotProtein(groundStateFold,hyperParams)

optimized_circuit = ProteinConfigAnsatz(angles);
sv = simulate(optimized_circuit);
histogram(sv,[1:8 10],Threshold=0.02)

reg = "us-east-1";
bucketPath = "s3://amazon-braket-mathworks/doc-examples";
device = quantum.backend.QuantumDeviceAWS("IonQ Device",S3Path=bucketPath,Region=reg)

task = run(optimized_circuit,device,NumShots=1000);
wait(task)

results = fetchOutput(task);
histogram(results,[1:8 10],Threshold=0.02)

