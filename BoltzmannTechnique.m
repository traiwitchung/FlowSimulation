%Boltzmann technique

load('geopack.mat');

A = 150;

L = 10e-6;

c = (L)/2; %half of voxel size

distgeo = L*bwdist(geopack,'euclidean');

for i = 1:A
    for j = 1:A
        for k = 1:A
            
if geopack(i,j,k) == 1
    D1(i,j,k) = 0;
    
else
    
    D1(i,j,k) = (distgeo(i,j,k) - c).^2 ;
    
end
        end
    end
end
    

clear geosmall
clear distgeo
clear i
clear j
clear k
clear c
clear L
clear A


%%
% D3 = D1(:);
% meandd = mean(mean(mean(D3(D3~=0))))
