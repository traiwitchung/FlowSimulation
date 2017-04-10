A = 300;

R = 1;

cc = (R)/2; %half of voxel size

distgeo = R*bwdist(papersample,'euclidean');
distgeo = round(distgeo/0.5)*0.5;

    

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

Dmax = zeros(13,18);

for i = 1:13
 for j = 1:18
  
  
  Dmax(i,j) = distgeo(i,j);
  
  end
    
 end
  




  for t = 1: 10
      
        for i = 2:13-1
        for j = 2:18-1
        
            
      if Dmax(i,j) ~= 0
          
      if distgeo(i,j) < distgeo(i+1,j);
          
          Dmax(i,j) = Dmax(i+1,j);
          
      end
      
      if distgeo(i,j) < distgeo(i-1,j);
          
              Dmax(i,j) = Dmax(i-1,j);
      end 
      
       if distgeo(i,j) < distgeo(i+1,j+1);
          
          Dmax(i,j) = Dmax(i+1,j+1);
          
       end
      
       if distgeo(i,j) < distgeo(i-1,j-1);
          
              Dmax(i,j) = Dmax(i-1,j-1);
      end 
          
      if distgeo(i,j) < distgeo(i,j+1);
          
              Dmax(i,j) = Dmax(i,j+1);
      end 
            
      if distgeo(i,j) < distgeo(i,j-1);
          
              Dmax(i,j) = Dmax(i,j-1);
      end
      
      if distgeo(i,j) < distgeo(i-1,j+1);
          
              Dmax(i,j) = Dmax(i-1,j+1);
      end 
      
      if distgeo(i,j) < distgeo(i+1,j-1);
          
              Dmax(i,j) = Dmax(i+1,j-1);
      end
              
    
                    end
                   end
                  end
        end
 
              
  %%
  for i = 1:size(distgeo,1)-1;
    
    figure(100)
    temp = Dmax(:,:);
    imagesc(squeeze(temp));
    axis equal
    colorbar
    drawnow
    
  end
  
  %%
% L = 300;
% mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
% G=cartGrid([L L L]);
% 
%     figure
%     plotCellData(G, distgeo(:));
%     s.EdgeColor = 'none';
%     colorbar;
%     view(3);
    
    drawnow
