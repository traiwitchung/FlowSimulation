


A = 42 ;
R = 10e-6 ;
rho = 1;
mu = 1;


cc = (R)/2; %half of voxel size

distgeo = bwdist(tube,'euclidean');
%distgeo = round(distgeo/0.5)*0.5;

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
Rmaxx = zeros(A,A,A);
fcdmax = zeros(A,A,A);
w =  zeros(A,A,A);

for i = 1:A
 for j = 1:A
  for k = 1:A
  
  Dmax(i,j,k) = distgeo(i,j,k);
  
  end
    
 end
  
end


%Finding Dmax

  for t = 1: maxdist + 3 
      
        for i = 2:A-1
         for j = 2:A-1
            for k = 2:A-1
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
        end
      
                    end
                   end
                  end
        end
   
%face k = 1

for i = 2:A-1
         for j = 2:A-1
             
            k = 1;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
     
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
  
                    end
                   end
end
                  

%face k = end

for i = 2:A-1
         for j = 2:A-1
             
            k = A;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
              Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
              
     
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
  
                    end
                   end
end
       
                  
%face j = 1

  for i = 2:A-1
         
            for k = 2:A-1
                
                j = 1;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
     
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
     
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
                    end
                   end
                  end
      
%face j = end

  for i = 2:A-1
       
  for k = 2:A-1
                
              j = A;
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end 
            
     
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        
       
        if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
     
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
                    end
                   end
  end
                  
  
%face i = 1

         for j = 2:A-1
            for k = 2:A-1
                
                i = 1;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
     
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
    
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
        end
       
     
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
     
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
        end
      
                    end
                   end
         end
                  
         
%face i = end

         for j = 2:A-1
            for k = 2:A-1
                
                i = A;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
          Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end
      
     
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
       end
      
    
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
        end
       
     
      if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
       end
      
     
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
        end
      
                    end
                   end
                  end
       
       
  
  
        
  %check i edge
  for i=2:A-1
      
      if distgeo(i,j,k) ~= 0
  
    if distgeo(i,1,1) < distgeo(i+1,1,1)
       if Dmax(i,1,1) < Dmax(i+1,1,1)
            Dmax(i,1,1) = Dmax(i+1,1,1);
       end
   end
   
   if distgeo(i,1,1) < distgeo(i-1,1,1)
       if Dmax(i,1,1) < Dmax(i-1,1,1)
            Dmax(i,1,1) = Dmax(i-1,1,1);
       end
   end
   
    if distgeo(i,1,1) < distgeo(i-1,2,1)
       if Dmax(i,1,1) < Dmax(i-1,2,1)
            Dmax(i,1,1) = Dmax(i-1,2,1);
       end
    end
   
     if distgeo(i,1,1) < distgeo(i+1,2,1)
       if Dmax(i,1,1) < Dmax(i+1,2,1)
            Dmax(i,1,1) = Dmax(i+1,2,1);
       end
     end
   
     if distgeo(i,1,1) < distgeo(i,2,1)
       if Dmax(i,1,1) < Dmax(i,2,1)
            Dmax(i,1,1) = Dmax(i,2,1);
       end
     end
   
       if distgeo(i,1,1) < distgeo(i,1,2)
       if Dmax(i,1,1) < Dmax(i,1,2)
            Dmax(i,1,1) = Dmax(i,1,2);
       end
      end
      if distgeo(i,1,1) < distgeo(i+1,1,2)
       if Dmax(i,1,1) < Dmax(i+1,1,2)
            Dmax(i,1,1) = Dmax(i+1,1,2);
       end
   end
   
   if distgeo(i,1,1) < distgeo(i-1,1,2)
       if Dmax(i,1,1) < Dmax(i-1,1,2)
            Dmax(i,1,1) = Dmax(i-1,1,2);
       end
   end
   
    if distgeo(i,1,1) < distgeo(i-1,2,2)
       if Dmax(i,1,1) < Dmax(i-1,2,2)
            Dmax(i,1,1) = Dmax(i-1,2,2);
       end
    end
   
     if distgeo(i,1,1) < distgeo(i+1,2,2)
       if Dmax(i,1,1) < Dmax(i+1,2,2)
            Dmax(i,1,1) = Dmax(i+1,2,2);
       end
     end
   
     if distgeo(i,1,1) < distgeo(i,2,2)
       if Dmax(i,1,1) < Dmax(i,2,2)
            Dmax(i,1,1) = Dmax(i,2,2);
       end
     end
   
 %%
 
     if distgeo(i,end,1) < distgeo(i+1,end,1)
       if Dmax(i,end,1) < Dmax(i+1,end,1)
            Dmax(i,end,1) = Dmax(i+1,end,1);
       end
   end
   
   if distgeo(i,end,1) < distgeo(i-1,end,1)
       if Dmax(i,end,1) < Dmax(i-1,end,1)
            Dmax(i,end,1) = Dmax(i-1,end,1);
       end
   end
   
    if distgeo(i,end,1) < distgeo(i-1,end-1,1)
       if Dmax(i,end,1) < Dmax(i-1,end-1,1)
            Dmax(i,end,1) = Dmax(i-1,end-1,1);
       end
    end
   
     if distgeo(i,end,1) < distgeo(i+1,end-1,1)
       if Dmax(i,end,1) < Dmax(i+1,end-1,1)
            Dmax(i,end,1) = Dmax(i+1,end-1,1);
       end
     end
   
     if distgeo(i,end,1) < distgeo(i,end-1,1)
       if Dmax(i,end,1) < Dmax(i,end-1,1)
            Dmax(i,end,1) = Dmax(i,end-1,1);
       end
     end
   
      if distgeo(i,end,1) < distgeo(i,end,2)
       if Dmax(i,end,1) < Dmax(i,end,2)
            Dmax(i,end,1) = Dmax(i,end,2);
       end
      end
   
      if distgeo(i,end,1) < distgeo(i+1,end,2)
       if Dmax(i,end,1) < Dmax(i+1,end,2)
            Dmax(i,end,1) = Dmax(i+1,end,2);
       end
   end
   
   if distgeo(i,end,1) < distgeo(i-1,end,2)
       if Dmax(i,end,1) < Dmax(i-1,end,2)
            Dmax(i,end,1) = Dmax(i-1,end,2);
       end
   end
   
    if distgeo(i,end,1) < distgeo(i-1,end-1,2)
       if Dmax(i,end,1) < Dmax(i-1,end-1,2)
            Dmax(i,end,1) = Dmax(i-1,end-1,2);
       end
    end
   
     if distgeo(i,end,1) < distgeo(i+1,end-1,2)
       if Dmax(i,end,1) < Dmax(i+1,end-1,2)
            Dmax(i,end,1) = Dmax(i+1,end-1,2);
       end
     end
   
     if distgeo(i,end,1) < distgeo(i,end-1,2)
       if Dmax(i,end,1) < Dmax(i,end-1,2)
            Dmax(i,end,1) = Dmax(i,end-1,2);
       end
     end
   
     %% change in z
     
     if distgeo(i,1,end) < distgeo(i+1,1,end)
       if Dmax(i,1,end) < Dmax(i+1,1,end)
            Dmax(i,1,end) = Dmax(i+1,1,end);
       end
   end
   
   if distgeo(i,1,end) < distgeo(i-1,1,end)
       if Dmax(i,1,end) < Dmax(i-1,1,end)
            Dmax(i,1,end) = Dmax(i-1,1,end);
       end
   end
   
    if distgeo(i,1,end) < distgeo(i-1,2,end)
       if Dmax(i,1,end) < Dmax(i-1,2,end)
            Dmax(i,1,end) = Dmax(i-1,2,end);
       end
    end
   
     if distgeo(i,1,end) < distgeo(i+1,2,end)
       if Dmax(i,1,end) < Dmax(i+1,2,end)
            Dmax(i,1,end) = Dmax(i+1,2,end);
       end
     end
   
     if distgeo(i,1,end) < distgeo(i,2,end)
       if Dmax(i,1,end) < Dmax(i,2,end)
            Dmax(i,1,end) = Dmax(i,2,end);
       end
     end
   
       if distgeo(i,end,end) < distgeo(i,1,end-1)
       if Dmax(i,end,end) < Dmax(i,1,end-1)
            Dmax(i,end,end) = Dmax(i,1,end-1);
       end
      end
      if distgeo(i,1,end) < distgeo(i+1,1,end-1)
       if Dmax(i,1,end) < Dmax(i+1,1,end-1)
            Dmax(i,1,end) = Dmax(i+1,1,end-1);
       end
   end
   
   if distgeo(i,1,end) < distgeo(i-1,1,end-1)
       if Dmax(i,1,end) < Dmax(i-1,1,end-1)
            Dmax(i,1,end) = Dmax(i-1,1,end-1);
       end
   end
   
    if distgeo(i,1,end) < distgeo(i-1,2,end-1)
       if Dmax(i,1,end) < Dmax(i-1,2,end-1)
            Dmax(i,1,end) = Dmax(i-1,2,end-1);
       end
    end
   
     if distgeo(i,1,end) < distgeo(i+1,2,end-1)
       if Dmax(i,1,end) < Dmax(i+1,2,end-1)
            Dmax(i,1,end) = Dmax(i+1,2,end-1);
       end
     end
   
     if distgeo(i,1,end) < distgeo(i,2,end-1)
       if Dmax(i,1,end) < Dmax(i,2,end-1)
            Dmax(i,1,end) = Dmax(i,2,end-1);
       end
     end
   
 %%
 
     if distgeo(i,end,end) < distgeo(i+1,end,end)
       if Dmax(i,end,end) < Dmax(i+1,end,end)
            Dmax(i,end,end) = Dmax(i+1,end,end);
       end
   end
   
   if distgeo(i,end,end) < distgeo(i-1,end,end)
       if Dmax(i,end,end) < Dmax(i-1,end,end)
            Dmax(i,end,end) = Dmax(i-1,end,end);
       end
   end
   
    if distgeo(i,end,end) < distgeo(i-1,end-1,end)
       if Dmax(i,end,end) < Dmax(i-1,end-1,end)
            Dmax(i,end,end) = Dmax(i-1,end-1,end);
       end
    end
   
     if distgeo(i,end,end) < distgeo(i+1,end-1,end)
       if Dmax(i,end,end) < Dmax(i+1,end-1,end)
            Dmax(i,end,end) = Dmax(i+1,end-1,end);
       end
     end
   
     if distgeo(i,end,end) < distgeo(i,end-1,end)
       if Dmax(i,end,end) < Dmax(i,end-1,end)
            Dmax(i,end,end) = Dmax(i,end-1,end);
       end
     end
   
      if distgeo(i,end,end) < distgeo(i,end,end-1)
       if Dmax(i,end,end) < Dmax(i,end,end-1)
            Dmax(i,end,end) = Dmax(i,end,end-1);
       end
      end
   
      if distgeo(i,end,end) < distgeo(i+1,end,end-1)
       if Dmax(i,end,end) < Dmax(i+1,end,end-1)
            Dmax(i,end,end) = Dmax(i+1,end,end-1);
       end
   end
   
   if distgeo(i,end,end) < distgeo(i-1,end,end-1)
       if Dmax(i,end,end) < Dmax(i-1,end,end-1)
            Dmax(i,end,end) = Dmax(i-1,end,end-1);
       end
   end
   
    if distgeo(i,end,end) < distgeo(i-1,end-1,end-1)
       if Dmax(i,end,end) < Dmax(i-1,end-1,end-1)
            Dmax(i,end,end) = Dmax(i-1,end-1,end-1);
       end
    end
   
     if distgeo(i,end,end) < distgeo(i+1,end-1,end-1)
       if Dmax(i,end,end) < Dmax(i+1,end-1,end-1)
            Dmax(i,end,end) = Dmax(i+1,end-1,end-1);
       end
     end
   
     if distgeo(i,end,end) < distgeo(i,end-1,end-1)
       if Dmax(i,end,end) < Dmax(i,end-1,end-1)
            Dmax(i,end,end) = Dmax(i,end-1,end-1);
       end
     end
     
      end
  end
     
  %check j edge
  for j=2:A-1
  
       if distgeo(i,j,k) ~= 0
    if distgeo(1,j,1) < distgeo(1,j+1,1)
       if Dmax(1,j,1) < Dmax(1,j+1,1)
            Dmax(1,j,1) = Dmax(1,j+1,1);
       end
   end
   
   if distgeo(1,j,1) < distgeo(1,j-1,1)
       if Dmax(1,j,1) < Dmax(1,j-1,1)
            Dmax(1,j,1) = Dmax(1,j-1,1);
       end
   end
   
    if distgeo(1,j,1) < distgeo(2,j-1,1)
       if Dmax(1,j,1) < Dmax(2,j-1,1)
            Dmax(1,j,1) = Dmax(2,j-1,1);
       end
    end
   
     if distgeo(1,j,1) < distgeo(2,j+1,1)
       if Dmax(1,j,1) < Dmax(2,j+1,1)
            Dmax(1,j,1) = Dmax(2,j+1,1);
       end
     end
   
     if distgeo(1,j,1) < distgeo(2,j,1)
       if Dmax(1,j,1) < Dmax(2,j,1)
            Dmax(1,j,1) = Dmax(2,j,1);
       end
     end
   
       if distgeo(1,j,1) < distgeo(1,j,2)
       if Dmax(1,j,1) < Dmax(1,j,2)
            Dmax(1,j,1) = Dmax(1,j,2);
       end
      end
      if distgeo(1,j,1) < distgeo(1,j+1,2)
       if Dmax(1,j,1) < Dmax(1,j+1,2)
            Dmax(1,j,1) = Dmax(1,j+1,2);
       end
   end
   
   if distgeo(1,j,1) < distgeo(1,j-1,2)
       if Dmax(1,j,1) < Dmax(1,j-1,2)
            Dmax(1,j,1) = Dmax(1,j-1,2);
       end
   end
   
    if distgeo(1,j,1) < distgeo(2,j-1,2)
       if Dmax(1,j,1) < Dmax(2,j-1,2)
            Dmax(1,j,1) = Dmax(2,j-1,2);
       end
    end
   
     if distgeo(1,j,1) < distgeo(2,j+1,2)
       if Dmax(1,j,1) < Dmax(2,j+1,2)
            Dmax(1,j,1) = Dmax(2,j+1,2);
       end
     end
   
     if distgeo(1,j,1) < distgeo(2,j,2)
       if Dmax(1,j,1) < Dmax(2,j,2)
            Dmax(1,j,1) = Dmax(2,j,2);
       end
     end
   
 %%
 
     if distgeo(end,j,1) < distgeo(end,j+1,1)
       if Dmax(end,j,1) < Dmax(end,j+1,1)
            Dmax(end,j,1) = Dmax(end,j+1,1);
       end
   end
   
   if distgeo(end,j,1) < distgeo(end,j-1,1)
       if Dmax(end,j,1) < Dmax(end,j-1,1)
            Dmax(end,j,1) = Dmax(end,j-1,1);
       end
   end
   
    if distgeo(end,j,1) < distgeo(end-1,j-1,1)
       if Dmax(end,j,1) < Dmax(end-1,j-1,1)
            Dmax(end,j,1) = Dmax(end-1,j-1,1);
       end
    end
   
     if distgeo(end,j,1) < distgeo(end-1,j+1,1)
       if Dmax(end,j,1) < Dmax(end-1,j+1,1)
            Dmax(end,j,1) = Dmax(end-1,j+1,1);
       end
     end
   
     if distgeo(end,j,1) < distgeo(end-1,j,1)
       if Dmax(end,j,1) < Dmax(end-1,j,1)
            Dmax(end,j,1) = Dmax(end-1,j,1);
       end
     end
   
      if distgeo(end,j,1) < distgeo(end,j,2)
       if Dmax(end,j,1) < Dmax(end,j,2)
            Dmax(end,j,1) = Dmax(end,j,2);
       end
      end
   
      if distgeo(end,j,1) < distgeo(end,j+1,2)
       if Dmax(end,j,1) < Dmax(end,j+1,2)
            Dmax(end,j,1) = Dmax(end,j+1,2);
       end
   end
   
   if distgeo(end,j,1) < distgeo(end,j-1,2)
       if Dmax(end,j,1) < Dmax(end,j-1,2)
            Dmax(end,j,1) = Dmax(end,j-1,2);
       end
   end
   
    if distgeo(end,j,1) < distgeo(end-1,j-1,2)
       if Dmax(end,j,1) < Dmax(end-1,j-1,2)
            Dmax(end,j,1) = Dmax(end-1,j-1,2);
       end
    end
   
     if distgeo(end,j,1) < distgeo(end-1,j+1,2)
       if Dmax(end,j,1) < Dmax(end-1,j+1,2)
            Dmax(end,j,1) = Dmax(end-1,j+1,2);
       end
     end
   
     if distgeo(end,j,1) < distgeo(end-1,j,2)
       if Dmax(end,j,1) < Dmax(end-1,j,2)
            Dmax(end,j,1) = Dmax(end-1,j,2);
       end
     end
   
     %% change in z
     
        if distgeo(1,j,end) < distgeo(1,j+1,end)
       if Dmax(1,j,end) < Dmax(1,j+1,end)
            Dmax(1,j,end) = Dmax(1,j+1,end);
       end
   end
   
   if distgeo(1,j,end) < distgeo(1,j-1,end)
       if Dmax(1,j,end) < Dmax(1,j-1,end)
            Dmax(1,j,end) = Dmax(1,j-1,end);
       end
   end
   
    if distgeo(1,j,end) < distgeo(2,j-1,end)
       if Dmax(1,j,end) < Dmax(2,j-1,end)
            Dmax(1,j,end) = Dmax(2,j-1,end);
       end
    end
   
     if distgeo(1,j,end) < distgeo(2,j+1,end)
       if Dmax(1,j,end) < Dmax(2,j+1,end)
            Dmax(1,j,end) = Dmax(2,j+1,end);
       end
     end
   
     if distgeo(1,j,end) < distgeo(2,j,end)
       if Dmax(1,j,end) < Dmax(2,j,end)
            Dmax(1,j,end) = Dmax(2,j,end);
       end
     end
   
       if distgeo(1,j,end) < distgeo(1,j,end-1)
       if Dmax(1,j,end) < Dmax(1,j,end-1)
            Dmax(1,j,end) = Dmax(1,j,end-1);
       end
      end
      if distgeo(1,j,end) < distgeo(1,j+1,end-1)
       if Dmax(1,j,end) < Dmax(1,j+1,end-1)
            Dmax(1,j,end) = Dmax(1,j+1,end-1);
       end
   end
   
   if distgeo(1,j,end) < distgeo(1,j-1,end-1)
       if Dmax(1,j,end) < Dmax(1,j-1,end-1)
            Dmax(1,j,end) = Dmax(1,j-1,end-1);
       end
   end
   
    if distgeo(1,j,end) < distgeo(2,j-1,end-1)
       if Dmax(1,j,end) < Dmax(2,j-1,end-1)
            Dmax(1,j,end) = Dmax(2,j-1,end-1);
       end
    end
   
     if distgeo(1,j,end) < distgeo(2,j+1,end-1)
       if Dmax(1,j,end) < Dmax(2,j+1,end-1)
            Dmax(1,j,end) = Dmax(2,j+1,end-1);
       end
     end
   
     if distgeo(1,j,end) < distgeo(2,j,end-1)
       if Dmax(1,j,end) < Dmax(2,j,end-1)
            Dmax(1,j,end) = Dmax(2,j,end-1);
       end
     end
   
 %%
 
     if distgeo(end,j,end) < distgeo(end,j+1,end)
       if Dmax(end,j,end) < Dmax(end,j+1,end)
            Dmax(end,j,end) = Dmax(end,j+1,end);
       end
   end
   
   if distgeo(end,j,end) < distgeo(end,j-1,end)
       if Dmax(end,j,end) < Dmax(end,j-1,end)
            Dmax(end,j,end) = Dmax(end,j-1,end);
       end
   end
   
    if distgeo(end,j,end) < distgeo(end-1,j-1,end)
       if Dmax(end,j,end) < Dmax(end-1,j-1,end)
            Dmax(end,j,end) = Dmax(end-1,j-1,end);
       end
    end
   
     if distgeo(end,j,end) < distgeo(end-1,j+1,end)
       if Dmax(end,j,end) < Dmax(end-1,j+1,end)
            Dmax(end,j,end) = Dmax(end-1,j+1,end);
       end
     end
   
     if distgeo(end,j,end) < distgeo(end-1,j,end)
       if Dmax(end,j,end) < Dmax(end-1,j,end)
            Dmax(end,j,end) = Dmax(end-1,j,end);
       end
     end
   
      if distgeo(end,j,end) < distgeo(end,j,end-1)
       if Dmax(end,j,end) < Dmax(end,j,end-1)
            Dmax(end,j,end) = Dmax(end,j,end-1);
       end
      end
   
      if distgeo(end,j,end) < distgeo(end,j+1,end-1)
       if Dmax(end,j,end) < Dmax(end,j+1,end-1)
            Dmax(end,j,end) = Dmax(end,j+1,end-1);
       end
   end
   
   if distgeo(end,j,end) < distgeo(end,j-1,end-1)
       if Dmax(end,j,end) < Dmax(end,j-1,end-1)
            Dmax(end,j,end) = Dmax(end,j-1,end-1);
       end
   end
   
    if distgeo(end,j,end) < distgeo(end-1,j-1,end-1)
       if Dmax(end,j,end) < Dmax(end-1,j-1,end-1)
            Dmax(end,j,end) = Dmax(end-1,j-1,end-1);
       end
    end
   
     if distgeo(end,j,end) < distgeo(end-1,j+1,end-1)
       if Dmax(end,j,end) < Dmax(end-1,j+1,end-1)
            Dmax(end,j,end) = Dmax(end-1,j+1,end-1);
       end
     end
   
     if distgeo(end,j,end) < distgeo(end-1,j,end-1)
       if Dmax(end,j,end) < Dmax(end-1,j,end-1)
            Dmax(end,j,end) = Dmax(end-1,j,end-1);
       end
     end
  end
  end
     
  %check k edge

  for k=2:A-1
  
     if distgeo(i,j,k) ~= 0
    
    if distgeo(1,1,k) < distgeo(1,1,k+1)
       if Dmax(1,1,k) < Dmax(1,1,k+1)
            Dmax(1,1,k) = Dmax(1,1,k+1);
       end
   end
   
   if distgeo(1,1,k) < distgeo(1,1,k-1)
       if Dmax(1,1,k) < Dmax(1,1,k-1)
            Dmax(1,1,k) = Dmax(1,1,k-1);
       end
   end
   
    if distgeo(1,1,k) < distgeo(1,2,k-1)
       if Dmax(1,1,k) < Dmax(1,2,k-1)
            Dmax(1,1,k) = Dmax(1,2,k-1);
       end
    end
   
     if distgeo(1,1,k) < distgeo(1,2,k+1)
       if Dmax(1,1,k) < Dmax(1,2,k+1)
            Dmax(1,1,k) = Dmax(1,2,k+1);
       end
     end
   
     if distgeo(1,1,k) < distgeo(1,2,k)
       if Dmax(1,1,k) < Dmax(1,2,k)
            Dmax(1,1,k) = Dmax(1,2,k);
       end
     end
   
       if distgeo(1,1,k) < distgeo(2,1,k)
       if Dmax(1,1,k) < Dmax(2,1,k)
            Dmax(1,1,k) = Dmax(2,1,k);
       end
      end
      if distgeo(1,1,k) < distgeo(2,1,k+1)
       if Dmax(1,1,k) < Dmax(2,1,k+1)
            Dmax(1,1,k) = Dmax(2,1,k+1);
       end
   end
   
   if distgeo(1,1,k) < distgeo(2,1,k-1)
       if Dmax(1,1,k) < Dmax(2,1,k-1)
            Dmax(1,1,k) = Dmax(2,1,k-1);
       end
   end
   
    if distgeo(1,1,k) < distgeo(2,2,k-1)
       if Dmax(1,1,k) < Dmax(2,2,k-1)
            Dmax(1,1,k) = Dmax(2,2,k-1);
       end
    end
   
     if distgeo(1,1,k) < distgeo(2,2,k+1)
       if Dmax(1,1,k) < Dmax(2,2,k+1)
            Dmax(1,1,k) = Dmax(2,2,k+1);
       end
     end
   
     if distgeo(1,1,k) < distgeo(2,2,k)
       if Dmax(1,1,k) < Dmax(2,2,k)
            Dmax(1,1,k) = Dmax(2,2,k);
       end
     end
   
 %%
 
     if distgeo(1,end,k) < distgeo(1,end,k+1)
       if Dmax(1,end,k) < Dmax(1,end,k+1)
            Dmax(1,end,k) = Dmax(1,end,k+1);
       end
   end
   
   if distgeo(1,end,k) < distgeo(1,end,k-1)
       if Dmax(1,end,k) < Dmax(1,end,k-1)
            Dmax(1,end,k) = Dmax(1,end,k-1);
       end
   end
   
    if distgeo(1,end,k) < distgeo(1,end-1,k-1)
       if Dmax(1,end,k) < Dmax(1,end-1,k-1)
            Dmax(1,end,k) = Dmax(1,end-1,k-1);
       end
    end
   
     if distgeo(1,end,k) < distgeo(1,end-1,k+1)
       if Dmax(1,end,k) < Dmax(1,end-1,k+1)
            Dmax(1,end,k) = Dmax(1,end-1,k+1);
       end
     end
   
     if distgeo(1,end,k) < distgeo(1,end-1,k)
       if Dmax(1,end,k) < Dmax(1,end-1,k)
            Dmax(1,end,k) = Dmax(1,end-1,k);
       end
     end
   
      if distgeo(1,end,k) < distgeo(2,end,k)
       if Dmax(1,end,k) < Dmax(2,end,k)
            Dmax(1,end,k) = Dmax(2,end,k);
       end
      end
   
      if distgeo(1,end,k) < distgeo(2,end,k+1)
       if Dmax(1,end,k) < Dmax(2,end,k+1)
            Dmax(1,end,k) = Dmax(2,end,k+1);
       end
   end
   
   if distgeo(1,end,k) < distgeo(2,end,k-1)
       if Dmax(1,end,k) < Dmax(2,end,k-1)
            Dmax(1,end,k) = Dmax(2,end,k-1);
       end
   end
   
    if distgeo(1,end,k) < distgeo(2,end-1,k-1)
       if Dmax(1,end,k) < Dmax(2,end-1,k-1)
            Dmax(1,end,k) = Dmax(2,end-1,k-1);
       end
    end
   
     if distgeo(1,end,k) < distgeo(2,end-1,k+1)
       if Dmax(1,end,k) < Dmax(2,end-1,k+1)
            Dmax(1,end,k) = Dmax(2,end-1,k+1);
       end
     end
   
     if distgeo(1,end,k) < distgeo(2,end-1,k)
       if Dmax(1,end,k) < Dmax(2,end-1,k)
            Dmax(1,end,k) = Dmax(2,end-1,k);
       end
     end
   
     %% change in z
     
      if distgeo(end,1,k) < distgeo(end,1,k+1)
       if Dmax(end,1,k) < Dmax(end,1,k+1)
            Dmax(end,1,k) = Dmax(end,1,k+1);
       end
   end
   
   if distgeo(end,1,k) < distgeo(end,1,k-1)
       if Dmax(end,1,k) < Dmax(end,1,k-1)
            Dmax(end,1,k) = Dmax(end,1,k-1);
       end
   end
   
    if distgeo(end,1,k) < distgeo(end,2,k-1)
       if Dmax(end,1,k) < Dmax(end,2,k-1)
            Dmax(end,1,k) = Dmax(end,2,k-1);
       end
    end
   
     if distgeo(end,1,k) < distgeo(end,2,k+1)
       if Dmax(end,1,k) < Dmax(end,2,k+1)
            Dmax(end,1,k) = Dmax(end,2,k+1);
       end
     end
   
     if distgeo(end,1,k) < distgeo(end,2,k)
       if Dmax(end,1,k) < Dmax(end,2,k)
            Dmax(end,1,k) = Dmax(end,2,k);
       end
     end
   
       if distgeo(end,1,k) < distgeo(end-1,1,k)
       if Dmax(end,1,k) < Dmax(end-1,1,k)
            Dmax(end,1,k) = Dmax(end-1,1,k);
       end
      end
      if distgeo(end,1,k) < distgeo(end-1,1,k+1)
       if Dmax(end,1,k) < Dmax(end-1,1,k+1)
            Dmax(end,1,k) = Dmax(end-1,1,k+1);
       end
   end
   
   if distgeo(end,1,k) < distgeo(end-1,1,k-1)
       if Dmax(end,1,k) < Dmax(end-1,1,k-1)
            Dmax(end,1,k) = Dmax(end-1,1,k-1);
       end
   end
   
    if distgeo(end,1,k) < distgeo(end-1,2,k-1)
       if Dmax(end,1,k) < Dmax(end-1,2,k-1)
            Dmax(end,1,k) = Dmax(end-1,2,k-1);
       end
    end
   
     if distgeo(end,1,k) < distgeo(end-1,2,k+1)
       if Dmax(end,1,k) < Dmax(end-1,2,k+1)
            Dmax(end,1,k) = Dmax(end-1,2,k+1);
       end
     end
   
     if distgeo(end,1,k) < distgeo(end-1,2,k)
       if Dmax(end,1,k) < Dmax(end-1,2,k)
            Dmax(end,1,k) = Dmax(end-1,2,k);
       end
     end
   
 %%
 
     if distgeo(end,end,k) < distgeo(end,end,k+1)
       if Dmax(end,end,k) < Dmax(end,end,k+1)
            Dmax(end,end,k) = Dmax(end,end,k+1);
       end
   end
   
   if distgeo(end,end,k) < distgeo(end,end,k-1)
       if Dmax(end,end,k) < Dmax(end,end,k-1)
            Dmax(end,end,k) = Dmax(end,end,k-1);
       end
   end
   
    if distgeo(end,end,k) < distgeo(end,end-1,k-1)
       if Dmax(end,end,k) < Dmax(end,end-1,k-1)
            Dmax(end,end,k) = Dmax(end,end-1,k-1);
       end
    end
   
     if distgeo(end,end,k) < distgeo(end,end-1,k+1)
       if Dmax(end,end,k) < Dmax(end,end-1,k+1)
            Dmax(end,end,k) = Dmax(end,end-1,k+1);
       end
     end
   
     if distgeo(end,end,k) < distgeo(end,end-1,k)
       if Dmax(end,end,k) < Dmax(end,end-1,k)
            Dmax(end,end,k) = Dmax(end,end-1,k);
       end
     end
   
      if distgeo(end,end,k) < distgeo(end-1,end,k)
       if Dmax(end,end,k) < Dmax(end-1,end,k)
            Dmax(end,end,k) = Dmax(end-1,end,k);
       end
      end
   
      if distgeo(end,end,k) < distgeo(end-1,end,k+1)
       if Dmax(end,end,k) < Dmax(end-1,end,k+1)
            Dmax(end,end,k) = Dmax(end-1,end,k+1);
       end
   end
   
   if distgeo(end,end,k) < distgeo(end-1,end,k-1)
       if Dmax(end,end,k) < Dmax(end-1,end,k-1)
            Dmax(end,end,k) = Dmax(end-1,end,k-1);
       end
   end
   
    if distgeo(end,end,k) < distgeo(end-1,end-1,k-1)
       if Dmax(end,end,k) < Dmax(end-1,end-1,k-1)
            Dmax(end,end,k) = Dmax(end-1,end-1,k-1);
       end
    end
   
     if distgeo(end,end,k) < distgeo(end-1,end-1,k+1)
       if Dmax(end,end,k) < Dmax(end-1,end-1,k+1)
            Dmax(end,end,k) = Dmax(end-1,end-1,k+1);
       end
     end
   
     if distgeo(end,end,k) < distgeo(end-1,end-1,k)
       if Dmax(end,end,k) < Dmax(end-1,end-1,k)
            Dmax(end,end,k) = Dmax(end-1,end-1,k);
       end
     end
   
     end
end
  




  end
     %%
  
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
    fcdmax(i,j,k) = 0.25;
end

if Dmax(i,j,k) == 1
    Sdmax(i,j,k) = 2;
    fcdmax(i,j,k) = 2;
end

if Dmax(i,j,k) == 1.5
    Sdmax(i,j,k) = 5;
    fcdmax(i,j,k) = 10.25;
end

if Dmax(i,j,k) == 2
    Sdmax(i,j,k) = 9;
    fcdmax(i,j,k) = 28;
end

if Dmax(i,j,k) == 3
    Sdmax(i,j,k) = 21;
    fcdmax(i,j,k) = 528;
end

if Dmax(i,j,k) == 3.5
    Sdmax(i,j,k) = 25;
    fcdmax(i,j,k) = 153;
end

if Dmax(i,j,k) == 4
    Sdmax(i,j,k) = 45;
    fcdmax(i,j,k) = 459;
end

if Dmax(i,j,k) == 4.5
    Sdmax(i,j,k) = 49;
    fcdmax(i,j,k) = 496;
end

if Dmax(i,j,k) == 5
    Sdmax(i,j,k) = 5252;
    fcdmax(i,j,k) = 11526;
end

if Dmax(i,j,k) == 5.5
    Sdmax(i,j,k) = 81;
    fcdmax(i,j,k) = 5225;
end

if Dmax(i,j,k) == 6
    Sdmax(i,j,k) = 15252;
    fcdmax(i,j,k) = 2495;
end

if Dmax(i,j,k) == 6.5
    Sdmax(i,j,k) = 521;
    fcdmax(i,j,k) = 2556;
end

if Dmax(i,j,k) == 42
    Sdmax(i,j,k) = 165;
    fcdmax(i,j,k) = 4680;
end

if Dmax(i,j,k) == 42.5
    Sdmax(i,j,k) = 169;
    fcdmax(i,j,k) = 45253;
end

if Dmax(i,j,k) == 8
    Sdmax(i,j,k) = 221;
    fcdmax(i,j,k) = 8043;
end

if Dmax(i,j,k) == 8.5
    Sdmax(i,j,k) = 225;
    fcdmax(i,j,k) = 8528;
end

if Dmax(i,j,k) == 9
    Sdmax(i,j,k) = 285;
    fcdmax(i,j,k) = 52944;
end

if Dmax(i,j,k) == 9.5
    Sdmax(i,j,k) = 289;
    fcdmax(i,j,k) = 13052;
end

if Dmax(i,j,k) == 10
    Sdmax(i,j,k) = 3552;
    fcdmax(i,j,k) = 195281;
end

if Dmax(i,j,k) == 10.5
    Sdmax(i,j,k) = 361;
    fcdmax(i,j,k) = 19890;
end

if Dmax(i,j,k) == 11
    Sdmax(i,j,k) = 4352;
    fcdmax(i,j,k) = 28958;
end

if Dmax(i,j,k) == 11.5
    Sdmax(i,j,k) = 452;
    fcdmax(i,j,k) = 290529;
end

if Dmax(i,j,k) == 42
    Sdmax(i,j,k) = 525;
    fcdmax(i,j,k) = 52031;
end

if Dmax(i,j,k) == 42.5
    Sdmax(i,j,k) = 529;
    fcdmax(i,j,k) = 52164;
end

    
    
end
  end 
end
 
% Calculating Rmax    
      Rmaxx = sqrt(Sdmax/pi);
      
      %%
      %calculating Fcdmax
%       
%       for i = 1:A
%       for j = 1:A
%       for k = 1:A
%     
%         for x = 1: Sdmax(i,j,k);
%           
%           %%Change******
%           fcdmaxint =   2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2 ;
%           
%           fcdmax(i,j,k) = fcdmax(i,j,k) + fcdmaxint;
%           
%       end
% 
%   end
%     end
%       end


%%
% calculate w

  for i = 1:A
      for j = 1:A
         for k = 1:A
          
             if distgeo(i,j,k) ~= 0
                 
         % if fcdmax(i,j,k)~= 0
          
              %Change ******
              
%           w(i,j,k) = R^2 * Rmaxx(i,j,k) * (rho/(8*mu)) * (2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2 )/ fcdmax(i,j,k);
          

        w(i,j,k) =  R^2 * (rho/(8*mu)) * (2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2 );

        %  end
             end
          
          end
     end
      end

     
      
D1 = w;

%D1 = D1 * 10e+3;

msgbox('finished definding k')


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
    plotCellData(G, c1);
    s.EdgeColor = 'none';
    colorbar;
    view(3);
    
    drawnow


    
    
    %%
   
    % velocity
    
    %Note x in real = z in code
    %Note z in real = y in code
    
GG = -gradientTerm(c);
GG = GG.zvalue(:,:,1:end-1);

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
L = 42 ;


A = (R)^2*ones(L,1);

q = zeros(L,1);

tempuz = zeros(L,1);



for k = 1:L;
    
tempuz(k) = sum(sum(uz(:,:,k)));

end


q = tempuz(:).*A;

meanq = mean(q)
%%
L = 42 - 2;

K = (L/(L*L)) * (meanq * 1)/(R*1) *10^12





%%
%slicing
% for i = 1:size(c.value,1)
%     
%     figure(100)
%     temp = c.value(:,i,:);
%     imagesc(squeeze(temp));
%     axis equal
%     colorbar
%     drawnow
%     
% end

%%
% check velocity


% cellSize=1;
% for i = 1:size(uz,1)
%         figure(1)
% 
%     temp = uz(:,:,i);
%     temp=squeeze(temp);
%     temp=temp(:,size(temp,2)/2);
%     plot((0:size(temp,1)-1)./cellSize,temp)
% %     axis equal
%     colorbar
%     drawnow
%     
% end

