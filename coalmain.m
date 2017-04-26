%%
%define k using BoltzmannTechnique

%Boltzmann technique

%a = msgbox('start')

%%

%only for tube
%load('tube.mat');


%%
A = 42;

R = 10e-6;

c = (R)/2; %half of voxel size

distgeo = R*bwdist(tube,'euclidean');

for i = 1:A
    for j = 1:A
        for k = 1:A
            
if tube(i,j,k) == 1
    D1(i,j,k) = 0;
    
else
    
    D1(i,j,k) = (distgeo(i,j,k) - c).^2 ;
    
end
        end
    end
end
    
%D1 = D1*10e+3;

clear geosmall
clear distgeo
clear i
clear j
clear k
clear c
clear L
clear A
 

clear tube
% D3 = D1(:);
% meandd = mean(mean(mean(D3(D3~=0))))
% 
% clear D3

%msgbox('finished definding k')

%% Diffusion equation heavily optimised for large grids
%clc; clear;

%FVToolStartUp()

%  load('tube')

% L = size(geosmall(:,1,1))

%% Define the domain and create a mesh structure
L = 42;  % domain length
Nx = L; % number of cells
m = createMesh3D(Nx,Nx,Nx,L,L,L);

%%

L = 42;
mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([L L L]);

%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
% BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
% BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
% BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=1; % top boundary
% BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary
BC.front.a(:) = 0; BC.front.b(:)=1; BC.front.c(:)=0; % front boundary
BC.back.a(:) = 0; BC.back.b(:)=1; BC.back.c(:)=1; % back boundary
%% define the transfer coeffs
% D1 = rand(m.dims);
% D1(D1<0.4)= 0;
% D1(D1>0)= 0.5;
% D1(D1>0.8)= 10*D1(D1>0.8);
% D1 = 1

% for i = 1:L
%     for j = 1:L
%         for k = 1:L
% if geosmall(i,j,k) == 1
%     D1(i,j,k) = 0;
%     
% else
%     
%     D1(i,j,k) = 1;
%     
% end
% 
%         end
%     end
% end
% 
%     


  %D2 = D1(:);

%BT produces D1
%MRST produces D2


D = zeros(42.^3,1);
D(G.cells.indexMap)= D1(:);
D = reshape(D,m.dims);%rand(m.dims);

clear G
% D = 1

D = createCellVariable(m, D);
D.value=ndSparse(D.value,size(D.value));


    clear D1
    clear D2

%%

    Dface = harmonicMean(D);
    Dface.xvalue(isnan(Dface.xvalue))=0;
    Dface.yvalue(isnan(Dface.yvalue))=0;
    Dface.zvalue(isnan(Dface.zvalue))=0;
    %%
    D = diffusionTerm(Dface);
    
%     %Source Term
%     SourceTerm = createCellVariable(m,ones(m.dims)) ;
%     F = constantSourceTerm(SourceTerm);
    
    [M, RHS] = boundaryCondition(BC);
    RHS = RHS; %+ F ;
    M = D+M;
    
    
    
    %%
 % clear all the no need data
    clear D
%     clear SourceTerm
%     clear F
   % clear geo2
    clear tube
    clear i
    clear j
    clear k
    clear D1
    clear D2
    clear BC
    clear Nx
   % clear L
%     M = gpuArray(M);  
%     RHS=gpuArray(RHS);
%     c = gmres(M,RHS,[],[],100000);




%%
%Solveing PDE
 %msgbox('start solving')

    c = solvePDE(m, M, RHS);
    
    
    clear RHS
    clear m
    clear M
    
%%    
    %c.value(c.value==inf) = 0;
    
 %visualizeCells(c)

    c1 = c.value(2:end-1,2:end-1,2:end-1);
    c1 = full(c1);
    c1=c1(:);
    
%     
%     c2 = c.value(2:end-1,10,10);
%     c2 = c2(:);
    
%% %    %


% startup
% clc; clear;
L = 42;
mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([L L L]);

% 
    figure
    plotCellData(G, c1(:));
    s.EdgeColor = 'none';
    colorbar;
    view(3);
    
    drawnow

    
    
    %%
    
    % velocity
    
    %Note x in real = z in code
    %Note z in real = y in code


    u = Dface.*-gradientTerm(c)/R;

    uzvalue1 = u.zvalue(:,:,1:end-1);
    uzvalue = full(uzvalue1);
    
    uz = uzvalue;
    uz(isnan(uz))= 0;
    

%     uz = uz(:);
%     figure
%     plotCellData(G, uz)
%     colorbar
%     view(3)
%     
%     drawnow

    %%
% find flow rate of each layer
L = 42;


A = (R)^2*ones(L,1);

q = zeros(L,1);

tempuz = zeros(L,1);



for k = 1:L;
    
tempuz(k) = sum(sum(uz(:,:,k)));

end


q = tempuz(:).*A;

meanq = mean(q);
%%

L = L 

K = (L/(L*L)) * (meanq * 1)/(R*1) *10^12

%Kxx=1/((42*R).^2).*sum(sum(sum(full(uz))))

% check1 = uz(:,:,3);
% check2 = sum(sum(check1));
% a = sum(uz(:,:,1))
