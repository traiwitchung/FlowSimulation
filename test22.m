

Lx = 102
Ly = 4
tube = zeros(Lx,Ly);
shape = 4;

% 
%     for y =1:L
%         for x = 1:L
%             
%             if (x- 21)^2 + (y- 21)^2 <= 20^2 ;
%                 tube(x,y,:) = 0 ;
%                 
%             else
%                 
%                 tube(x,y,:) = 1;
%                 
%             end
%         end
%     end

%%

% 
    for y =1:Ly
        for x = 1:Lx
            
           tube(:,1) = 1;
           
           tube(:,1,:) = 1;
           
           tube(:,Ly) = 1;
           
          tube(:,Ly,:) = 1;
          
           
        end
    end
    
    clear x
    clear y
%     
%     mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
%     G=cartGrid([L L]);
% %
%     figure
%     plotCellData(G, tube(:));
%     

 %%
 R = 50e-6;
 
% 
distgeo = bwdist(tube,'euclidean');
maxdist = max(distgeo(:));


distgeo2 = R*bwdist(tube,'euclidean');
%distgeo = bwdist(tube,'euclidean');
%distgeo = round(distgeo/0.5)*0.5;



    

%% 
%visualise

% for i = 1:size(distgeo,1)
%     
%     figure(100)
%     temp = distgeo(i,:,:);
%     imagesc(squeeze(temp));
%     axis equal
%     colorbar
%     drawnow
%     
% end

%%
%define initial Dmax

Dmax = zeros(Lx,Ly);

for i = 1:Lx
 for j = 1:Ly
  
  
  Dmax(i,j) = distgeo(i,j);
  
  end
    
 end
  




  for t = 1: round(maxdist) +1
      
        for i = 2:Lx-1
        for j = 2:Ly-1
        
            
      if Dmax(i,j) ~= 0
          
      if distgeo(i,j) < distgeo(i+1,j);
          if Dmax(i,j) < Dmax(i+1,j);
          Dmax(i,j) = Dmax(i+1,j);
          end
          
      end
      
      if distgeo(i,j) < distgeo(i-1,j);
          if Dmax(i,j) < Dmax(i-1,j);
              Dmax(i,j) = Dmax(i-1,j);
          end
      end 
      
       if distgeo(i,j) < distgeo(i+1,j+1);
          if Dmax(i,j) < Dmax(i+1,j+1);
          Dmax(i,j) = Dmax(i+1,j+1);
          end
          
       end
      
       if distgeo(i,j) < distgeo(i-1,j-1);
          if Dmax(i,j) < Dmax(i-1,j-1);
              Dmax(i,j) = Dmax(i-1,j-1);
          end
      end 
          
      if distgeo(i,j) < distgeo(i,j+1);
          if Dmax(i,j) < Dmax(i,j+1);
              Dmax(i,j) = Dmax(i,j+1);
          end
      end 
            
      if distgeo(i,j) < distgeo(i,j-1);
          if Dmax(i,j) < Dmax(i,j-1);
              Dmax(i,j) = Dmax(i,j-1);
          end
      end
      
      if distgeo(i,j) < distgeo(i-1,j+1);
          if Dmax(i,j) < Dmax(i-1,j+1);
              Dmax(i,j) = Dmax(i-1,j+1);
          end
      end 
      
      if distgeo(i,j) < distgeo(i+1,j-1);
          if Dmax(i,j) < Dmax(i+1,j-1);
              Dmax(i,j) = Dmax(i+1,j-1);
          end
      end
              
    
                    end
                   end
                  end
  end
        
  Dmax(1,:) = Dmax(2,:);
  Dmax(Lx,:) = Dmax(Lx-1,:);
  
    %%
    %%
% calculate w
rho = 1;
mu = 1;

  for i = 1:Lx
      for j = 1:Ly
        
          
             if distgeo(i,j) ~= 0
                 
         
        w(i,j) = shape * R^2 * (rho/(8*mu)) * (2 * Dmax(i,j) * distgeo(i,j) - (distgeo(i,j))^2 );

       
             else
                 
                 w(i,j) = 0;
                 
             end
          
          end
     end
      

     
      
D1 = w;

    

    %%
    
mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([Lx Ly]);
    
Nx = Lx; % number of cells
m = createMesh2D(Lx,Ly, Lx,Ly);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure

BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
% BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=0; % top boundary
% BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary

%%
D = zeros(Lx,Ly);
D(G.cells.indexMap)= D1(:);
D = reshape(D,m.dims);

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
    
%       mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
%       G=cartGrid([Lx Ly]);
%         figure(10)
%         subplot(2,1,1)
%     plotCellData(G, c1);
%     
%     subplot(2,1,2)
%     plotCellData(G, tube(:));
    
    %%
    Dfacex = full(Dface.xvalue);
    
    gradient = -gradientTerm(c)
    
    u = Dface.*-gradientTerm(c)/1e-6;  %special case
    
    ux = u.xvalue(1:end-1,:);
    uy = u.yvalue(:,1:end-1);
    
    uxx = full(u.xvalue);
    ux = full(ux);
    ux(isnan(ux)) = 0;
    uy = full(uy);
    
%     figure(11)
%     plotCellData(G, ux(:));
   % figure(12)
   % plotCellData(G, uy(:));
    
    %%
    
    Uprofile = ux(Lx/2,:);
    Uprofile = Uprofile(:);
    
    
    %%
    
% A = (R)^2*ones(L,1);
% q = zeros(L,1);
% tempuz = zeros(L,1);
%     
% for k = 1:L;
%     
% tempuz(k) = sum(sum(ux(k,:)));
% 
% end
% 
% q = tempuz(:).*A;
% 
% meanq = mean(q)
% 
% K = (L/(L*L)) * (meanq * 1)/(R*1) 
    
    