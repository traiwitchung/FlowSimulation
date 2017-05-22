%profile -memory on;
%instruction
%put sample into workspace
%change sample, A and R
%done ez!


R = 10e-6/1 ;
rho = 1;
mu = 1;
Shape = 1;
%shape 1 is tube, 4 is fractures
sample = geopack;

A = size(sample,1);

poro = sum(sample(:)==0)/size(sample,1)^3

distgeo = bwdist(sample,'euclidean');
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

w =  zeros(A,A,A);

for i = 1:A
 for j = 1:A
  for k = 1:A
  
  Dmax(i,j,k) = distgeo(i,j,k);
  
  end
    
 end
  
end


%Finding Dmax

  for t = 1: round(maxdist)  
      
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
  %checking Dmax
  
  
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
% calculate w

  for i = 1:A
      for j = 1:A
         for k = 1:A
          
             if distgeo(i,j,k) ~= 0
                 

        w(i,j,k) =  Shape * R^2 * (rho/(8*mu)) * (2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2 );

             end
          
          end
     end
      end

     
      
D1 = w;

disp('finished definding k')



%% Define the domain and create a mesh structure

L = A;  % domain length
Nx = A; % number of cells
m = createMesh3D(Nx,Nx,Nx,L,L,L);

%%

L = A;
mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([L L L]);

%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure

%x direction
% BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
% BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary

%y direction
% BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=1; % top boundary
% BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary

%z direction
BC.front.a(:) = 0; BC.front.b(:)=1; BC.front.c(:)=0; % front boundary
BC.back.a(:) = 0; BC.back.b(:)=1; BC.back.c(:)=1; % back boundary
%% define the transfer coeffs
%put local conductivity into the cells

D = zeros(A.^3,1);
D(G.cells.indexMap)= D1(:);
D = reshape(D,m.dims);

clear G


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
    
    %preparing matrices
    
    D = diffusionTerm(Dface);
    
    [M, RHS] = boundaryCondition(BC);
    
    M = D+M;
    
    
    
    %%
 % clear all the no need data
    clear D
%     clear SourceTerm
%     clear F
   % clear geo2
   % clear sample
    clear i
    clear j
    clear k
    clear D1
    clear D2
    clear BC
    clear Nx
   % clear L




%%
%Solveing PDE

 
 
disp('start solving')

   c = solvePDE(m, M, RHS);
    
 % c = thomas(M,RHS);
    
    clear RHS
    clear m
    clear M
    
%%    
    %c.value(c.value==inf) = 0;
    
 %visualizeCells(c)

% c = reshape(c,152,152,152);
 
    c1 = c.value(2:end-1,2:end-1,2:end-1);
    c1 = full(c1);
    c1=c1(:);
    
%     

    
%% 


% startup

%plotting pressure profile


mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([A A A]);

    figure
    plotCellData(G, c1);
    s.EdgeColor = 'none';
    colorbar;
    view(3);
    
    drawnow


    
    
    %%
   
    disp('calculating K')
    
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
L = A ;


Area = (R)^2*ones(L,1);

q = zeros(L,1);

tempuz = zeros(L,1);



for k = 1:L;
    
tempuz(k) = sum(sum(uz(:,:,k)));

end


q = tempuz(:).*Area;

meanq = mean(q)
%%
L = A ;

K = (L/(L*L))*(meanq*1)/(R*1)*10^12 %Darcy





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

