function ansatz = ProteinConfigAnsatz(parameters)
% Create the circuit ansatz for a 7 amino acid neuropeptide (10 qubit circuit).
parameters = reshape(parameters, [2 9]);
%hGate([1:7 9 10])
gates = [
         ryGate([1:7 9 10], parameters(1,:))
         cxGate(1:4, 2:5)
         cxGate(5,10)
         cxGate(10,9)
         cxGate(9,8)
         cxGate(8,9)
         cxGate(9,8)
         cxGate([8 7 6], [7 6 1])
         ryGate([1:8 10], parameters(2,:))
        ];

ansatz = quantumCircuit(gates);
end




