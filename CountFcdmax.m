%create cube

clc; clear;


L = 23;
tube = zeros(L,L);



    for y =1:L
        for x = 1:L
            
           tube(1,:) = 1;
           
           tube(:,1) = 1;
           
           tube(L,:) = 1;
           
           tube(:,L) = 1;
           
           
           
   
           
           
        end
    end
% 
%     tube(2,2) = 1;
%     tube(L-1,L-1) = 1;
%     tube(2,L-1) = 1;
%     tube(L-1,2) = 1;

    
g2 = tube(:);


mrstModule add incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
G=cartGrid([L L]);
% 
    figure(1)
    plotCellData(G, g2);
    s.EdgeColor = 'none';
    colorbar;
    %view(3);
    
    drawnow
    
    clear G
    clear g2
    clear x
    clear y
    clear z
     clear s
     clear L
     
     %%
     
     distgeo = bwdist(tube,'euclidean');
distgeo = round(distgeo/0.5)*0.5;

maxdist = max(distgeo(:));

%%

fcdmax = (2*maxdist*0.5 - 0.5^2)*(sum(distgeo(:) == 0.5)) ...
         + (2*maxdist*1.0 - 1.0^2)*(sum(distgeo(:) == 1.0)) ...
         + (2*maxdist*1.5 - 1.5^2)*(sum(distgeo(:) == 1.5)) ...
         + (2*maxdist*2.0 - 2.0^2)*(sum(distgeo(:) == 2.0)) ...
         + (2*maxdist*2.5 - 2.5^2)*(sum(distgeo(:) == 2.5)) ...
         + (2*maxdist*3.0 - 3.0^2)*(sum(distgeo(:) == 3.0)) ...
         + (2*maxdist*3.5 - 3.5^2)*(sum(distgeo(:) == 3.5)) ...
         + (2*maxdist*4.0 - 4.0^2)*(sum(distgeo(:) == 4.0)) ...
         + (2*maxdist*4.5 - 4.5^2)*(sum(distgeo(:) == 4.5)) ...
         + (2*maxdist*5.0 - 5.0^2)*(sum(distgeo(:) == 5.0)) ...
         + (2*maxdist*5.5 - 5.5^2)*(sum(distgeo(:) == 5.5)) ...
         + (2*maxdist*6.0 - 6.0^2)*(sum(distgeo(:) == 6.0)) ...
         + (2*maxdist*6.5 - 6.5^2)*(sum(distgeo(:) == 6.5)) ...
         + (2*maxdist*7.0 - 7.0^2)*(sum(distgeo(:) == 7.0)) ...
         + (2*maxdist*7.5 - 7.5^2)*(sum(distgeo(:) == 7.5)) ...
         + (2*maxdist*8.0 - 8.0^2)*(sum(distgeo(:) == 8.0)) ...
         + (2*maxdist*8.5 - 8.5^2)*(sum(distgeo(:) == 8.5)) ...
         + (2*maxdist*9.0 - 9.0^2)*(sum(distgeo(:) == 9.0)) ...
         + (2*maxdist*9.5 - 9.5^2)*(sum(distgeo(:) == 9.5)) ...
         + (2*maxdist*10.0 - 10.5^2)*(sum(distgeo(:) == 10.0)) ...
         + (2*maxdist*10.5 - 10.5^2)*(sum(distgeo(:) == 10.5)) ...
         + (2*maxdist*11.0 - 11.0^2)*(sum(distgeo(:) == 11.0)) ...
         + (2*maxdist*11.5 - 11.5^2)*(sum(distgeo(:) == 11.5)) ...
         + (2*maxdist*12.0 - 12.0^2)*(sum(distgeo(:) == 12.0)) ...
         + (2*maxdist*12.5 - 12.5^2)*(sum(distgeo(:) == 12.5));
         
     
    

     
     