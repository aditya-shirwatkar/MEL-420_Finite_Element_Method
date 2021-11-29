clc
clear
%Materials
E = 200e9; %Elasticity modulus (Pa)
v = 0.25; %Poisson coefficient 
t = 0.036; %Plate thickness (m)
%Planar stress
D = (E/(1-v^2))*[1 v 0;
                 v 1 0;
                 0 0 (1-v)/2];
%Planar strain    
% D = (E/((1+v)*(1-2*v)))*[1-v v 0;
%                          v 1-v 0;
%                          0 0 (1-v)];

l = 2; % (m)
b = 1; % (m)
Ne_x = 6;
Ne_y = 3;
% Ne = Ne_x + Ne_y;
% mesh = cell(Ne,1);
% for ex=1:Ne_x
%     for ey=0:Ne_y
%         mesh{ex+ey} = [(ex-1)*l/Ne_x,     ey*b/Ne_y;
%                        (ex-1)*l/Ne_x, (ey+1)*b/Ne_y;
%                        (ex)*l/Ne_x,   (ey+1)*b/Ne_y];
%         
%         plot(mesh{ex+ey}(:,1), mesh{ex+ey}(:,2), 'o');
%         hold on
%     end
% end

% y0yn = reshape(linspace(0,b,Nn_y), [Nn_y, 1]);
% x0xn = reshape(linspace(0,l,Nn_x), [Nn_x, 1]);
y0yn = linspace(0,b,Ne_y+1);
x0xn = linspace(0,l,Ne_x+1);
[X, Y] = meshgrid(x0xn, y0yn);
T = delaunayTriangulation([X(:) Y(:)]);
triplot(T)
hold on
for i=1:size(T.Points,1)
    text(T.Points(i,1),T.Points(i,2),int2str(i), 'Color', 'red')
    hold on
end
% model = createpde('structural','modal-planestress');
% geometryFromEdges(model,@lshapeg);