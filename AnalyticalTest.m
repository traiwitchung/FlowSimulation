%% Diffusion equation heavily optimised for large grids
%clc; clear;

FVToolStartUp()
%% Define the domain and create a mesh structure
L = 50;  % domain length
Nx = L; % number of cells
m = createMesh3D(Nx,Nx,Nx,L,L,L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
%BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=0; % top boundary
%BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary
%% define the transfer coeffs
% D1 = rand(m.dims);
% D1(D1<0.4)= 0;
% D1(D1>0.8)= 10*D1(D1>0.8);

    D1 = 1;
    D2 = D1(:);



D = createCellVariable(m, D1);
D.value=ndSparse(D.value,size(D.value));
%%

    D = harmonicMean(D);
    D.xvalue(isnan(D.xvalue))=0;
    D.yvalue(isnan(D.yvalue))=0;
    D.zvalue(isnan(D.zvalue))=0;
    %%
    D = diffusionTerm(D);
    [M, RHS] = boundaryCondition(BC);
    M = D+M;
    clear D
%     M = gpuArray(M);  
%     RHS=gpuArray(RHS);
%     c = gmres(M,RHS,[],[],100000);
%%
    c = solvePDE(m,M, RHS);
    
    
    
    %visualizeCells(c)

    
    c1 = c.value(2:end-1,2:end-1,2:end-1);
    c1=c1(:);
    
    
    c2 = c.value(2:end-1,L,L);
    c2 = c2(:);
    
%     
%     figure
%     plotCellData(G, c1)
%     colorbar
%     
%     drawnow
%     
    %%
    x = m.cellcenters.x;
    for x = 1:L
         analytical(x) = (-1/50)*(x-0.5) + 1;
         
    end
    
 x = m.cellcenters.x; %[1:50];
 
    analytical = analytical(:)
    
    
    figure
    plot(x, c2, 'o', x, analytical);
xlabel('X [unit]'); ylabel('P');
legend('Numerical', 'Analytical');


%%
% k = X + 1
%%%%%%%%%%%%

L = 50;  % domain length
Nx = L; % number of cells
m = createMesh3D(Nx,Nx,Nx,L,L,L);
% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
    
    %Define perm (D)
    D1 = ones(L, L, L);
    
    for a = 1:L
        
    D1(a,:,:) = a+0.5;
    
    end
    
    %D1 = D1(:);
    D2 = D1(:);



D = createCellVariable(m, D1);
D.value=ndSparse(D.value,size(D.value));
%%

    D = harmonicMean(D);
    D.xvalue(isnan(D.xvalue))=0;
    D.yvalue(isnan(D.yvalue))=0;
    D.zvalue(isnan(D.zvalue))=0;
    %%
    D = diffusionTerm(D);
    [M, RHS] = boundaryCondition(BC);
    M = D+M;
    clear D
%     M = gpuArray(M);  
%     RHS=gpuArray(RHS);
%     c = gmres(M,RHS,[],[],100000);
%%
    c = solvePDE(m,M, RHS);
    
    
    
    %visualizeCells(c)

    
    c1 = c.value(2:end-1,2:end-1,2:end-1);
    c1=c1(:);
    
    
    c2 = c.value(2:end-1,L,L);
    c2 = c2(:);
    
    
    %%
    
    
     
    for x = 1:L
         analytical(x) = (-1/log(51))*(log(x + 0.5)) + 1;
         
    end
    
 x = m.cellcenters.x %[1:50];
 
    analytical = analytical(:)
    
    
    figure
    plot(x, c2, 'o', x, analytical);
xlabel('X [unit]'); ylabel('P');
legend('Numerical', 'Analytical');
    

