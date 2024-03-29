clc
close all
clear all

%% for iteration start
for h = 0.9
    for t = 0.5:0.05:1
%% Generation of arc mesh
% Define the geometry of the arc
% h = 0.0;            % [mm] The distance between the center point(0,0) and the bottom of the arc center
% t = 0.9;            % [mm] The thickness of the arc
L = 15;             % [mm] The distance between (0,0) and left or right end point of the arc
R = h/2+L^2/(2*h);  % [mm] The radius of the arc
z = 5;              % [mm] The depth of the arc

% Define the mesh space for the arc
nx = 101;                              % The numder of nodes in the x-direction
ny = 3;                             % The number of nodes in the y-direction
nz = 11;                              % The number of nodes in the z-direction
num_nodes = nx*ny*nz;                % Total number of nodes
num_elems = (nx-1)*(ny-1)*(nz-1);    % Total number of elements
plot_mesh = 'no';

% Define the coordinate of the nodes for arc
X = linspace(-L,L,nx)';                            % Meshing in x-direction (X-coordinates of nodes)
y_t = zeros(nx,1);                                 % The top boundary of the arc +t/2
y_b = zeros(nx,1);                                 % The bottom boundary of the arc -t/2
Y = zeros(nx,ny);
y0 = zeros(nx,1);
Z = linspace(0,z,nz)';                             % Meshing in z-direction (Z-coordinates of nodes)

for i=1:nx
    if h > 0
        y0(i) = (L^2/(2*h)-h/2)-sqrt(R^2-X(i)^2);% Getting y-coordinates of arc's center line
    else
        y0(i) = 0;
    end
    y_t(i) = y0(i)+t/2;                          % Getting y-coordinates of arc's top line
    y_b(i) = y0(i)-t/2;                          % Getting y-coordinates of arc's bottom line
    Y(i,:) = linspace(y_b(i),y_t(i),ny);        % Meshing in y-direction (Y-coordinates of nodes)
end

% Generate the node coordinate for arc
Nodes = zeros(num_nodes, 3);     % Defining the coordinates of nodes in one matrix

for j1 = 1:nz
    for j2 = 1:ny
        for j3 = 1:nx
            Nodes((j2-1)*nx+j3+(j1-1)*nx*ny,:) = [X(j3),Y(j3,j2),Z(j1)];
        end
    end
end

%% Genertation of column mesh
% Define the mesh space for the column
nx_c = 4;
ny_c = 25;
nz_c = 6;
num_nodes_c = nx_c*ny_c*nz_c;
num_elems_c = (nx_c-1)*(ny_c-1)*(nz_c-1);

% Define the coordinate of the nodes for column
X_c = linspace(-1,1,nx_c)';
if h > 0
    y_t_c = (L^2/(2*h)-h/2)-sqrt(R^2-1)+t/2;
else
    y_t_c = t/2;
end
Y_c = linspace(y_t_c,18.25,ny_c)';
Z_c = linspace(0,z,nz_c)';

% Generate the node coordinate for arc
Nodes_c = zeros(num_nodes_c, 3);

for j1 = 1:nz_c
    for j2 = 1:ny_c
        for j3 = 1:nx_c
            Nodes_c((j2-1)*nx_c+j3+(j1-1)*nx_c*ny_c,:) = [X_c(j3),Y_c(j2),Z_c(j1)];
        end
    end
end

%% Generate the elements connected by nodes (8-nodes 3D element)
elem_connect = zeros(num_elems, 8);
elem_connect_c = zeros(num_elems_c, 8);

for k1 = 1:(nz-1)
    for k2 = 1:(ny-1)
        for k3 = 1:(nx-1)
            n1 = (k2-1)*nx+k3+(k1-1)*nx*ny;
            n2 = n1 + 1;
            n3 = n2 + nx;
            n4 = n3 - 1;
            n5 = n1 + nx*ny;
            n6 = n5 + 1;
            n7 = n6 + nx;
            n8 = n7 - 1;
            elem_ID = k3+(k2-1)*(nx-1)+(k1-1)*(nx-1)*(ny-1);
            elem_connect(elem_ID,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
        end
    end
end

for k1 = 1:(nz_c-1)
    for k2 = 1:(ny_c-1)
        for k3 = 1:(nx_c-1)
            n1 = (k2-1)*nx_c+k3+(k1-1)*nx_c*ny_c + num_nodes;
            n2 = n1 + 1;
            n3 = n2 + nx_c;
            n4 = n3 - 1;
            n5 = n1 + nx_c*ny_c;
            n6 = n5 + 1;
            n7 = n6 + nx_c;
            n8 = n7 - 1;
            elem_ID_c = k3+(k2-1)*(nx_c-1)+(k1-1)*(nx_c-1)*(ny_c-1);
            elem_connect_c(elem_ID_c,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
        end
    end
end

%% export the nodes and elements informations (Peparing for inserting files)
node_ids = (1:num_nodes)'; % Node IDs
elem_ids = (1:num_elems)'; % Element IDs

node_ids_c = (num_nodes+1:num_nodes+num_nodes_c)'; % Node IDs of column
elem_ids_c = (num_elems+1:num_elems+num_elems_c)'; % Element IDs of column

node_name = sprintf('nodes_insert_arc.txt');
elem_name = sprintf('elems_insert_arc.txt');

node_name_c = sprintf('nodes_insert_column.txt');
elem_name_c = sprintf('elems_insert_column.txt');

writematrix([node_ids, Nodes], node_name);
writematrix([elem_ids, elem_connect], elem_name);

writematrix([node_ids_c, Nodes_c], node_name_c);
writematrix([elem_ids_c, elem_connect_c], elem_name_c);

%% Plot the arc mesh
if strcmpi(plot_mesh,'yes')==1
  for jj1=1:num_elems/(nz-1)
        n1 = elem_connect(jj1,1);
        n2 = elem_connect(jj1,2);
        n3 = elem_connect(jj1,3);
        n4 = elem_connect(jj1,4);
        mesh_x = [Nodes(n1,1) Nodes(n2,1) Nodes(n3,1) Nodes(n4,1) Nodes(n1,1)];
        mesh_y = [Nodes(n1,2) Nodes(n2,2) Nodes(n3,2) Nodes(n4,2) Nodes(n1,2)];
        plot(mesh_x,mesh_y,'k','LineWidth',0.25)
        hold on;
  end
  axis equal;

  for jj2=1:num_elems_c/(nz_c-1)
        n1 = elem_connect_c(jj2,1) - num_nodes;
        n2 = elem_connect_c(jj2,2) - num_nodes;
        n3 = elem_connect_c(jj2,3) - num_nodes;
        n4 = elem_connect_c(jj2,4) - num_nodes;
        mesh_x_c = [Nodes_c(n1,1) Nodes_c(n2,1) Nodes_c(n3,1) Nodes_c(n4,1) Nodes_c(n1,1)];
        mesh_y_c = [Nodes_c(n1,2) Nodes_c(n2,2) Nodes_c(n3,2) Nodes_c(n4,2) Nodes_c(n1,2)];
        plot(mesh_x_c,mesh_y_c,'r','LineWidth',0.25)
        hold on;
  end

end

%% Write inp.files
% read the original .inp file
original_text = split(string(fileread('OR-AND-GATE-allin.inp')), newline);

% read the node information and element information
insert1 = split(string(fileread(node_name)), newline);
insert2 = split(string(fileread(elem_name)), newline);

insert3 = split(string(fileread(node_name_c)), newline);
insert4 = split(string(fileread(elem_name_c)), newline);

insert_text1 = insert1(1:num_nodes);
insert_text2 = insert2(1:num_elems);

insert_text3 = insert3(1:num_nodes_c);
insert_text4 = insert4(1:num_elems_c);

% insert the node and element information into the designated row
insert_index = [65009 65010 65019 65020]; %[*Node(arc) *Element(arc) *Node(column) *Element(column)]
insert_text = {insert_text1, insert_text2, insert_text3, insert_text4};
for ii = 1:length(insert_index)
    idx = insert_index(ii);
    original_text = [original_text(1:idx); insert_text{ii}; original_text(idx+1:end)];
    % update the row numbers after inserting new contents
    insert_index = insert_index + numel(insert_text{ii});
end

% output the modified .inp file
Gate_name = sprintf('OR-AND-GATE-h%de-02-t%de-02.inp',h*100,t*100);
fileID = fopen(Gate_name, 'w');
fprintf(fileID, '%s\n', original_text{:});
fclose(fileID);

%% For iteration end
    end
end