clc
clear
disp('-------------FEM Rod Solver---------------')


nodes = input('Enter the number of nodes:');
    if isempty(nodes)
        nodes = 5;
    end
ele = nodes-1;
E = input("Young's Modulus in (GPa): ") * 10^9;
    if isempty(E)
        E = 200 * 10^9;
    end
l = input('Length in (m):');
    if isempty(l)
        l = 3;
    end
den = input('Enter the density of material in kg/m^3:');
    if isempty(den)
        den = 7750;
    end
% x = zeros(nodes,1);
% %x is coordinate matrix
% for i = 2:nodes
%     x(i,1)=x(i-1,1)+l/ele;
% end
b = input('Enter 1 for uniform area else 0:');
    if isempty(b)
        b = 0;
    end
if(b==1)
    A1 = input('Area in mm^2:');
        if isempty(A1)
            A1 = 100*100 * 10^-6;
        end
    A2 = A1;
    mini = A1;
else
    A1 = input('Area at node 1 in mm^2:');
        if isempty(A1)
            A1 = 100*100 * 10^-6;
        end
    A2 = input('Area at last node in mm^2:');
        if isempty(A2)
            A2 = 50*50 * 10^-6;
        end
    mini = min(A1,A2);
end
step = (abs(A1-A2)/ele);
AM = linspace(A1, A2, nodes)'; % Area vector

K=zeros(nodes); % Global Stiffness Matrix
C=[1 -1;-1 1];
i=1;
while(i<=ele)
    Ke = ((AM(i,1)+AM(i+1,1))*E*ele/(2*l))*C;
    K(i:i+1,i:i+1) =  K(i:i+1,i:i+1)+Ke;
    i=i+1;
end
F0 = input('Enter the force matrix in kN')*1000;
    if isempty(F0)
        F0 = [5 1*10^3];
    end
F=zeros(nodes,1);
force_index = F0(:,1);
F(force_index) = F0(:,2);

U = zeros(nodes,1); % displacement
%b = input('Enter 1 for no boundary conditions for displacement at last node else 0');
%if b==0
 %   U(nodes,1)= input('Boundary Condition at last node :')
%end

constrains=input('\nEnter constrains:');
    if isempty(constrains)
        constrains = [1];
        solveable = setxor(1:nodes, constrains);
    end

U(solveable,1) = K(solveable, solveable) \ F(solveable,1);
strain = zeros(ele,1);
stress = zeros(ele,1);
for i=1:ele
    strain(i,1) = U(i+1,1)-U(i,1);
    stress(i,1) = strain(i,1)*E; 
end
%FE = zeros(nodes,1)
%for i=1:nodes-1
  % f = [(AM(i,1)*x(i,1)/2)+(AM(i,1)*x(i+1,1)/6)+(AM(i+1,1)*x(i,1)/6)+(AM(i+1,1)*x(i+1,1)/6);
 %       (AM(i,1)*x(i,1)/6)+(AM(i,1)*x(i+1,1)/6)+(AM(i+1,1)*x(i,1)/6)+(AM(i+1,1)*x(i+1,1)/2)]*2
   %f = f*(den*l^2/(2*ele))
%   FE(i:i+1,1) = FE(i:i+1,1)+f
%end
%FE = FE*(den*l^2)/48
X = zeros(nodes,1);
Y = zeros(nodes,1);

for i = 2:nodes
    X(i) = X(i-1)+l/nodes;
end

X_D = X+U;
Y_D = Y;
fprintf('---------Global Stiffness Matrix------------\n')
fprintf('EA(min)/l * \n')
disp(K/(E/(l)))
fprintf('---------Displacement in Nodes(mm)----------\n')
for i = 1:nodes
    fprintf('Node %d     %d\r\n',i,U(i,1))
end

fprintf('---------Strain in Elements----------\n')
for i = 1:ele
    fprintf('Element %d     %d\r\n',i,strain(i,1))
end

fprintf('---------Stress in Elements(N/mm^2)----------\n')
for i = 1:ele
    fprintf('Element %d     %d\r\n',i,stress(i,1))
end

patch(X_D,Y_D,U,'EdgeColor','interp','linewidth',5);
title('Resultant Displacement Plot')
colormap jet
cb = colorbar;