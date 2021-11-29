% Program for 2D beam analysis. The beam elements are 2 node_coordd with
% Hermitian Shape functions
clear
clc
fprintf('\n----------Program for 2D beam analysis----------\n')  

nn=input('\nEnter number of nodes:'); % 10
L=input('\nEnter length of the beam (m):'); % 3
node_coord = linspace(0,L,nn);
ndof = 2*nn; % Every node_coord has 2 dofs.
se = 4; % Size of elemental matrix, since 2 node_coordd beam element
ne = nn-1; % Number of elements
le = L/nn; % Length of the respective element

% Dofs which are constrained
constrains=input('\nEnter constrains:'); % constrains = [1,2,7,19,20];
solveable = setxor(1:ndof, constrains);

Ftype = input('\nEnter Ftype (udl or con):'); % 'udl'
    if Ftype == 'udl'
        % [element_1 element_2 F1; element_3 element_4 F2; ...]
        Felem = input('\nEnter Fudl:'); % [1 3 -800*10^3]
    elseif Ftype == 'con'
        % [x1 F1; x2 F2] : Note 0 < x < L
        Felem = input('\nEnter Fcon (Note 0 < x < L):'); % [1.2 -800*10^3; 0.8 -800*10^3]
    else
        Ftype = 'udl';
        Felem = [1 3 -800*10^3];
    end

b=input('\nEnter width of the beam (mm):')*10^-3; % 100
h=input('\nEnter height of the beam (mm):')*10^-3; % 100
I = b*h^3 / 12;    
E=input('\nEnter modulus of elasticity of the element (GPa):') * 10^9; % 210

if Ftype == 'udl'
    Fudl_elem = zeros(ne,1);
    Fudl_elem(Felem(:,1):Felem(:,2))= Felem(:,3);
elseif Ftype == 'con'
    Fcon_elem = zeros(ne,2);
    Felem_index = [];
    for i=1:size(Felem,1)
        Felem_index = [Felem_index; find(node_coord > Felem(i,1), 1, 'first')];
    end
    a0_star = (node_coord(Felem_index-1) + node_coord(Felem_index))/2;
    a1_star = (node_coord(Felem_index) - node_coord(Felem_index-1))/2;
    zeta = (Felem(:,1) - a0_star')./a1_star';
    Fcon_elem(Felem_index, 1)= zeta;
    Fcon_elem(Felem_index, 2)= Felem(:,2);
end

K = zeros(ndof,ndof); % Define the global stiffness matrix
w = zeros(ndof,1); % Define the displacement matrix
F = zeros(ndof,1); % Define the global load vector

% Elemental stiffness matrix
ke = (E*I/le^3)*[12, 6*le, -12, 6*le;
                 6*le, 4*le^2, -6*le, 2*le^2;
                 -12, -6*le, 12, -6*le;
                 6*le, 2*le^2, -6*le, 4*le^2];

% Loop retrive global stiffness matrix and load vector
for e=1:ne
    
    % Elemental load vector
    if Ftype == 'udl'
        P0 = Fudl_elem(e);
        fe = [P0*le/2;
              P0*(le^2)/12;
              P0*le/2;
              -P0*(le^2)/12;];
    elseif Ftype == 'con'
        zeta_s = Fcon_elem(e,1); 
        P = Fcon_elem(e,2);
        H_s = [1/4 * (2 - 3*zeta_s + zeta_s^3);
               le/8 * (1 - zeta_s - zeta_s^2 + zeta_s^3);
               1/4 * (2 + 3*zeta_s - zeta_s^3);
               le/8 * (-1 - zeta_s + zeta_s^2 + zeta_s^3)];
        fe = P*le/2 * H_s;
    end
    
    % Add ke and fe to the global matrix and load vector 
    conn_elem_global= [2*e-1, 2*e, 2*(e+1)-1, 2*(e+1)];
    K(conn_elem_global, conn_elem_global) = K(conn_elem_global, conn_elem_global) + ke;
    F(conn_elem_global, :) = F(conn_elem_global, :) + fe;
end

% Solve for displacement and slope
w(solveable) = K(solveable, solveable) \ F(solveable);
w_dw_data = [];
fprintf('\n----------Nodal Displacements and Slope---------\n');
for i=1:nn
    w_dw_data = [w_dw_data; i, w(2*i-1), w(2*i)];
end
disp(table(w_dw_data(:,1), w_dw_data(:,2), w_dw_data(:,3), ...
            'VariableNames', {'Node Number','w', 'dw/dx'}))

% Reactions at nodes
f=K*w-F;
reactions_data = [];
fprintf('\n----------Reactions----------\n');
for i=1:nn
    if ismember(2*i-1, constrains) || ismember(2*i, constrains)
        reactions_data = [reactions_data; i, f(2*i-1)/1000, f(2*i)/1000];
    end
end
disp(table(reactions_data(:,1), reactions_data(:,2), reactions_data(:,3), ...
            'VariableNames', {'Node Number', 'Force(kN)', 'Moment(kN-m)'}))

% Elemental bending moment & bending stress
bending_stress_moment_shear_force = [];
for e=1:ne
    we1 = w(2*e-1);
    dwe1 = w(2*e);
    we2 = w(2*(e+1)-1);
    dwe2 = w(2*(e+1));
    Ve = (8*E*I/le^3)*(1.5*we1 + 0.75*dwe1 - 1.5*we2 + 0.75*dwe2);
    elemental_data = [];
    % Bending Moment: Me = 4EI/le^2 diff(We, z, 2)
    for ze=-1:1
        d2We_dz2 = 3/2*ze* we1 + le/4*(-1 + 3*ze)*dwe1 + ...
                   -3/2*ze*we2 + le/4*(1 + 3*ze)*dwe2;
        Me =(4*E*I/le^2)*d2We_dz2;
        stress_max = Me*(h/2)/I;
        a0 = (node_coord(e) + node_coord(e+1))/2;
        a1 = (-node_coord(e) + node_coord(e+1))/2;
        dom = a0 + a1*ze;
        elemental_data = [elemental_data; dom, Me/1000, stress_max/10^6, Ve/1000];
    end
    bending_stress_moment_shear_force = [bending_stress_moment_shear_force; 
                                         elemental_data];
end

% Stress bending graph
N = size(bending_stress_moment_shear_force, 1);

figure(1)
subplot(3,1,1);
a = plot(bending_stress_moment_shear_force(:,1), bending_stress_moment_shear_force(:,2),'-r');
legend('Bending Moment (kN-m)')
subplot(3,1,2);
b = plot(bending_stress_moment_shear_force(:,1), bending_stress_moment_shear_force(:,3),'-g');
legend('Max Bending Stress (MPa)')
subplot(3,1,3);
c = plot(bending_stress_moment_shear_force(:,1), bending_stress_moment_shear_force(:,4),'-b');
legend('Shear Force (kN)')

% Graph
figure(2)
plot(node_coord,zeros(nn,1), 'black', node_coord, w(1:2:end), 'r-o');
hold on
for i=1:nn
    if ismember(2*i-1, constrains) || ismember(2*i, constrains)
        plot(node_coord(i), 0,... 
                's',...
                'LineWidth',2,...
                'MarkerSize',10,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor',[0.5,0.5,0.5])
    end
    hold on
end

for i=1:ne
    if Ftype == 'udl'
        if abs(Fudl_elem(i)) > 0
            text(node_coord(i), w(2*i-1),'\downarrow')
            text(0.75*node_coord(i)+0.25*node_coord(i+1), 0.75*w(2*i-1)+0.25*w(2*(i+1)-1),'\downarrow')
            text(0.5*node_coord(i)+0.5*node_coord(i+1), 0.5*w(2*i-1)+0.5*w(2*(i+1)-1), strcat('\downarrow ', num2str(Fudl_elem(i)/1000),' kN'))
            text(0.25*node_coord(i)+0.75*node_coord(i+1), 0.25*w(2*i-1)+0.75*w(2*(i+1)-1),'\downarrow')
            text(node_coord(i+1), w(2*(i+1)-1),'\downarrow')
        end
        hold on
    elseif Ftype == 'con'
        if abs(Fcon_elem(i,2)) > 0
            z = Fcon_elem(i,1);
            x = (node_coord(i)+node_coord(i+1))/2 + z*(node_coord(i+1)-node_coord(i))/2;
            h1 = 1/4 * (2 - 3*z + z^3);
            h2 = le/8 * (1 - z - z^2 + z^3);
            h3 = 1/4 * (2 + 3*z - z^3);
            h4 = le/8 * (-1 - z + z^2 + z^3);
            we = h1*w(2*i-1) + h3*w(2*(i+1)-1) + h2*w(2*i) + h4*w(2*(i+1));
            text(x, we, strcat('\downarrow ', num2str(Fcon_elem(i,2)/1000),' kN'))
        end
        hold on
    end
end
xlabel('x (m)')
ylabel('bending deflection (m)')
ylim([-max(abs(w(1:2:end))) max(abs(w(1:2:end)))])