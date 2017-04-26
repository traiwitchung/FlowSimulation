L = 102;

tube = zeros(L,L);



    for y =1:L
        for x = 1:L
            
           tube(:,1) = 1;
           
           %tube(:,1,:) = 1;
           
           tube(:,L) = 1;
           
          % tube(:,L,:) = 1;
          
           
        end
    end
    
    clear x
    clear y
    
%     mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
%     G=cartGrid([L L]);
% %
%     figure
%     plotCellData(G, tube(:));
%     






%%
R = 4e-6;

c = (R)/2; %half of voxel size

distgeo = R*bwdist(tube,'euclidean');

for i = 1:L
    for j = 1:L
                    
if tube(i,j) == 1
    D1(i,j) = 0;
    
else
    
    D1(i,j) = (distgeo(i,j) - c).^2 ;
    
end
        end
end

    %%
    
     
mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([L L]);
    
Nx = L; % number of cells
m = createMesh2D(Nx,Nx, L,L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure

BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
% BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=0; % top boundary
% BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary

%%
D = zeros(L.^2,1);
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
    Dface.xvalue(isnan(Dface.xvalue))= 0;
    Dface.yvalue(isnan(Dface.yvalue))= 0;
    
        D = diffusionTerm(Dface);
    
%     %Source Term
%     SourceTerm = createCellVariable(m,ones(m.dims)) ;
%     F = constantSourceTerm(SourceTerm);
    
    [M, RHS] = boundaryCondition(BC);
    %RHS = RHS; %+ F ;
    M = D+M;
    
       %%
 % clear all the no need data
    clear D
%     clear SourceTerm
%     clear F
   % clear geo2
   % clear tube
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
    
    c1 = c.value(2:end-1,2:end-1);
    c1 = full(c1);
    c1=c1(:);
    
    %%
    
      mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
      G=cartGrid([L L]);
        figure(10)
        subplot(2,1,1)
    plotCellData(G, c1);
    
    subplot(2,1,2)
    plotCellData(G, tube(:));
    
    %%
    Dfacex = full(Dface.xvalue);
    
    gradient = -gradientTerm(c)
    
    u = Dface.*-gradientTerm(c)/R;
    
    ux = u.xvalue(1:end-1,:);
    uy = u.yvalue(:,1:end-1);
    
    uxx = full(u.xvalue);
    ux = full(ux);
    ux(isnan(ux)) = 0;
    uy = full(uy);
    
    figure(11)
    plotCellData(G, ux(:));
    figure(12)
    plotCellData(G, uy(:));
    
    %%
    
    Uprofile = ux(L/2,:);
    Uprofile = Uprofile(:);
    
    
    %%
    
A = (R)^2*ones(L,1);
q = zeros(L,1);
tempuz = zeros(L,1);
    
for k = 1:L;
    
tempuz(k) = sum(sum(ux(k,:)));

end

q = tempuz(:).*A;

meanq = mean(q)

K = (L/(L*L)) * (meanq * 1)/(R*1) 