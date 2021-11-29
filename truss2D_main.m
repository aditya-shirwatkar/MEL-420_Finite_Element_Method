clear 
clc
disp('------% Program for 2D truss analysis %------')

% Give the co-ordinates of the nodes
node=input('\nEnter node:'); % [0 0; 1 0; 1 1; 2 0; 2 1; 3 0];

% Defines the node connection of the trusses
conn=input('\nEnter conn:'); % [1 2; 1 3; 2 3; 2 4; 4 5; 3 5; 3 4; 4 6; 5 6];
    
nn=size(node,1);
ndof=2*nn; % Every node has 2 dofs.
ne=size(conn,1); % Size of elemental matrix

% Dofs which are constrained
constrains=input('\nEnter constrains:'); % [1,2,9,10]
solveable = setxor(1:ndof, constrains);

% Give the external loading (kN) [node_number_1 F1x F1y; node_number_2 F2x F2y; ...]
Fnodes = input('\nEnter Fnodes:'); % [3 80*10^3 -300*10^3; 5 400*10^3 -400*10^3]

A=input('\nEnter cross-sectional area of the element (mm^2):') * 10^-6; % 200*10^-6
E=input('\nEnter modulus of elasticity of the element (GPa):') * 10^9; % 210.0 * 10^9;

F=zeros(2*nn,1);
force_index_x = 2*Fnodes(:,1)-1;
force_index_y = 2*Fnodes(:,1);
F(force_index_x)= Fnodes(:,2);
F(force_index_y)= Fnodes(:,3);

K = zeros(ndof,ndof); % Define the global stiffness matrix 
d = zeros(ndof,1); % Define the displacement matrix
% Loop retrive global stiffness matrix for 2D truss
for e=1:ne
    n1=conn(e,1); % 1st node for the respective element
    n2=conn(e,2); % 2nd node for the respective element
    
    x1=node(n1,1); % X1 co-ordinate of the 1st node
    y1=node(n1,2); % Y1 co-ordinate of the 1st node
    x2=node(n2,1); % X2 co-ordinate of the 2nd node
    y2=node(n2,2); % Y2 co-ordinate of 2nd node
    
    L=sqrt(((x2-x1))^2+((y2-y1)^2)); % Length of the respective truss
    
    l=(x2-x1)/L; % cos(theta);
    m=(y2-y1)/L; % sin(theta);
    
    % Elemental stiffness matrix
    ke = (E*A/L)*[l^2, m*l, -l^2, -m*l;
                  m*l, m^2, -m*l, -m^2;
                  -l^2, -m*l, l^2, m*l;
                  -m*l, -m^2, m*l, m^2];
    
    % Add ke to the global matrix according to the respective node
    % and respective displacement
    elemental2global= [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    K(elemental2global, elemental2global) = K(elemental2global, elemental2global) + ke;
end

% Solve the matrices to get the displacements at each node
d(solveable) = K(solveable, solveable) \ F(solveable);

displacements_data = [];
fprintf('\n----------Nodal Displacements (mm)---------\n');
for i=1:nn
    displacements_data = [displacements_data; i, d(2*i-1)*1000, d(2*i)*1000];
end
disp(table(displacements_data(:,1), displacements_data(:,2), displacements_data(:,3), ...
            'VariableNames', {'Node Number','X-Direction', 'Y-Direction'}))

% Reactions at nodes
f=K*d-F;
reactions_data = [];
fprintf('\n----------Reactions (KN)----------\n');
for i=1:nn
    if ismember(2*i-1, constrains) || ismember(2*i, constrains)
        reactions_data = [reactions_data; i, f(2*i-1)/1000, f(2*i)/1000];
    end
end
disp(table(reactions_data(:,1), reactions_data(:,2), reactions_data(:,3), ...
            'VariableNames', {'Node Number', 'X-Direction', 'Y-Direction'}))

% Calculate the strains and stresses
fprintf('\n----------Elemental strain & stress----------\n');
stress_strain_data = [];
for e=1:ne
    n1=conn(e,1);
    n2=conn(e,2);
    
    x1=node(n1,1);
    y1=node(n1,2);
    x2=node(n2,1);
    y2=node(n2,2);
    
    L=sqrt(((x2-x1))^2+((y2-y1)^2));
    l=(x2-x1)/L;    
    m=(y2-y1)/L;
    
    B=(1/L)*[-l, -m, l, m];
    
    elemental2global = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    strain = B * d(elemental2global);
    stress = E * strain;
    
    stress_strain_data = [stress_strain_data; e, strain*1000, stress/10^6];
end
disp(table(stress_strain_data(:,1), stress_strain_data(:,2), stress_strain_data(:,3), ...
            'VariableNames', {'Member Number', 'Strain (mm)', 'Stress (MPa)'}))