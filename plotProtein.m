function plotProtein(bitstring,hyperParams)
% Plot protein structure from the bitstring

% The input bitstring is expected to be of length 9, with the first 7 bits
% specifying turns in direction in the structure, and the last 2 bits
% specifying interactions between beads.

% Number of beads
N = length(hyperParams.protein);

% Construct 3D coordinates representing the 4 corners of a tetrahedron
% centered in 0. These represent the 4 directions each new bead might
% be added on in a tetrahedral grid.
turn2bead = ones(4,3);
turn2bead(2:4,:) = -1+2*eye(3);

% Construct complete bitstring by inserting input bitstring into
% the complete string mask.
completeBitstring = hyperParams.turn2qubit;
completeBitstring(completeBitstring=='q') = bitstring(1:hyperParams.numQubitsConfig);

 % Each pair of bits is converted into a number between 0 and 3
turns = bin2dec(reshape(completeBitstring,2,[])');

% Change direction in which we follow the tetrahedral grid in each bead
% on the line.
signs = (-1).^(0:N-1)';

% Compute placements of each bead.
beads = cumsum(signs.*[zeros(1,3);turn2bead(turns+1,:)]);

% Plot the beads connecting by lines, and add a text label to each
figure
plot3(beads(:,1),beads(:,2),beads(:,3),'.-','LineWidth',2,'MarkerSize',80,'SeriesIndex',1)
axis off

viewDir = -[6 4 1];
view(viewDir);

% Plot text, ensuring it is in front of the lines
beadsText = beads + 0.03*viewDir;
text(beadsText(:,1),beadsText(:,2),beadsText(:,3),hyperParams.protein', ...
    'FontWeight', 'bold', 'HorizontalAlignment','center')

% Whether there are interactions between any pair of beads is determined
% by the additional bits in the string.
interactions = [];
currInteractionQubit = hyperParams.numQubitsConfig+1;

for i=1:(N-5)
    for j=(i+5):2:N
        if bitstring(currInteractionQubit) == '1'
            interactions = [interactions;beads(i,:);beads(j,:);nan*ones(1,3)]; %#ok<AGROW>
        end
        currInteractionQubit = currInteractionQubit+1;
    end
end

if ~isempty(interactions)
    hold on
    plot3(interactions(:,1),interactions(:,2),interactions(:,3),'k--','LineWidth',2)
    hold off
    legend([hyperParams.protein+" Protein Structure";"Interactions"], "Location","southoutside")
else
    legend(hyperParams.protein+" Protein Structure", "Location","southoutside")
end
end

