


A = 150;
R = 10e-6 ;
rho = 1;
mu = 1;


cc = (R)/2; %half of voxel size

distgeo = bwdist(geopack,'euclidean');
distgeo = round(distgeo/0.5)*0.5;

maxdist = max(distgeo(:));
    

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

Dmax = zeros(A,A,A);
Sdmax = zeros(A,A,A);
Rmax = zeros(A,A,A);
fcdmax = zeros(A,A,A);
w =  zeros(A,A,A);

for i = 1:A
 for j = 1:A
  for k = 1:A
  
  Dmax(i,j,k) = distgeo(i,j,k);
  
  end
    
 end
  
end



  for t = 1: maxdist + 2;
      
        for i = 2:A-1
        for j = 2:A-1
        for k = 2:A-1
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k);
          
          Dmax(i,j,k) = Dmax(i+1,j,k);
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k);
          
              Dmax(i,j,k) = Dmax(i-1,j,k);
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k);
          
              Dmax(i,j,k) = Dmax(i,j+1,k);
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k);
          
              Dmax(i,j,k) = Dmax(i,j-1,k);
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1);
          
              Dmax(i,j,k) = Dmax(i,j,k+1);
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1);
          
             Dmax(i,j,k) = Dmax(i,j,k-1);
      end
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k+1);
          
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1);
          
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1);
          
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1);
          
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1);
          
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k+1);
          
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1);
          
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k+1);
          
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k);
          
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k);
          
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k);
          
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k);
          
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k-1);
          
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1);
          
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1);
          
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1);
          
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1);
          
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k-1);
          
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1);
          
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k-1);
          
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
        end
      
                    end
                   end
                  end
        end
  end
              
  
%   for i = 1:size(distgeo,1)-1;
%     
%     figure(100)
%     temp = Dmax(i,:,:);
%     imagesc(squeeze(temp));
%     axis equal
%     colorbar
%     drawnow
%     
%   end

%%
%Sdmax
  
for i = 1:A
 for j = 1:A
  for k = 1:A
  
if Dmax(i,j,k) == 0.5
    Sdmax(i,j,k) = 1;
end

if Dmax(i,j,k) == 1
    Sdmax(i,j,k) = 2;
end

if Dmax(i,j,k) == 1.5
    Sdmax(i,j,k) = 5;
end

if Dmax(i,j,k) == 2
    Sdmax(i,j,k) = 9;
end

if Dmax(i,j,k) == 3
    Sdmax(i,j,k) = 21;
end

if Dmax(i,j,k) == 3.5
    Sdmax(i,j,k) = 25;
end

if Dmax(i,j,k) == 4
    Sdmax(i,j,k) = 45;
end

if Dmax(i,j,k) == 4.5
    Sdmax(i,j,k) = 49;
end

if Dmax(i,j,k) == 5
    Sdmax(i,j,k) = 77;
end

if Dmax(i,j,k) == 5.5
    Sdmax(i,j,k) = 81;
end

if Dmax(i,j,k) == 6
    Sdmax(i,j,k) = 177;
end

if Dmax(i,j,k) == 6.5
    Sdmax(i,j,k) = 121;
end

if Dmax(i,j,k) == 7
    Sdmax(i,j,k) = 165;
end

if Dmax(i,j,k) == 7.5
    Sdmax(i,j,k) = 169;
end

if Dmax(i,j,k) == 8
    Sdmax(i,j,k) = 221;
end

if Dmax(i,j,k) == 8.5
    Sdmax(i,j,k) = 225;
end

if Dmax(i,j,k) == 9
    Sdmax(i,j,k) = 285;
end

if Dmax(i,j,k) == 9.5
    Sdmax(i,j,k) = 289;
end

if Dmax(i,j,k) == 10
    Sdmax(i,j,k) = 357;
end

if Dmax(i,j,k) == 10.5
    Sdmax(i,j,k) = 361;
end

if Dmax(i,j,k) == 11
    Sdmax(i,j,k) = 437;
end

if Dmax(i,j,k) == 11.5
    Sdmax(i,j,k) = 441;
end

if Dmax(i,j,k) == 12
    Sdmax(i,j,k) = 525;
end

if Dmax(i,j,k) == 12.5
    Sdmax(i,j,k) = 529;
end

    
    
end
  end 
end
 
% Calculating Rmax    
      Rmax = sqrt(Sdmax/pi);
      
      %%
      %calculating Fcdmax
      
      for i = 1:A
      for j = 1:A
      for k = 1:A
    
      for x = 1: Sdmax(i,j,k)
          
          fcdmaxint = 2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2;
          
          fcdmax(i,j,k) = fcdmax(i,j,k) + fcdmaxint;
          
      end

  end
    end
      end


%%
% calculate w

  for i = 1:A
      for j = 1:A
      for k = 1:A
          
          if fcdmax(i,j,k)~= 0
          
          w(i,j,k) = R^2 * Rmax(i,j,k) * (rho/(8*mu)) * (2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2)/ fcdmax(i,j,k);
          
          end
          
          
          
          end
     end
      end

     
      
D1 = w;


msgbox('finished definding k')


%% Define the domain and create a mesh structure
L = 150;  % domain length
Nx = L; % number of cells
m = createMesh3D(Nx,Nx,Nx,L,L,L);

%%

L = 150;
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


D = zeros(150.^3,1);
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
    clear geopack
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
 msgbox('start solving')

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
L = 150;
mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([L L L]);

%     figure
%     plotCellData(G, c1);
%     s.EdgeColor = 'none';
%     colorbar;
%     view(3);
%     
%     drawnow

    
    
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
L = 150;


A = (R)^2*ones(L,1);

q = zeros(L,1);

tempuz = zeros(L,1);



for k = 1:L;
    
tempuz(k) = sum(sum(uz(:,:,k)));

end


q = tempuz(:).*A;

meanq = mean(q)
%%


K = (L/(L*L)) * (meanq * 1)/(R*1) 

%Kxx=1/((150*R).^2).*sum(sum(sum(full(uz))))

% check1 = uz(:,:,3);
% check2 = sum(sum(check1));
% a = sum(uz(:,:,1))
