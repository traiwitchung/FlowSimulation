

%calc error between MRST and FVTool

% p = importdata('p.dat');
% p2 = p(:);

resSol.pressure = p2;


c1 = c1(G.cells.indexMap);

c1(isnan(c1)) = 0;



resSol.pressure(isnan(resSol.pressure)) = 0;

s = sum((c1 - resSol.pressure).^2)/size(c1,1)


%%

error = zeros(size(c1));

for i=1:size(c1,1)
  if resSol.pressure(i) ~= 0 && c1(i) ~= 0;
     error(i) = abs(c1(i)-resSol.pressure(i))./resSol.pressure(i);
  end
end
 

meanerror = mean(error)


%%
%Relative percentage different
for i=1:size(c1,1)
  if resSol.pressure(i) ~= 0 && c1(i) ~= 0;
     error(i) = 2*(c1(i)-resSol.pressure(i))./(abs(c1(i))+abs(resSol.pressure(i)));
  end
end
 


RPD = mean(error)


%%

% 
% 
%  figure
%     plot(x, c2, 'o', x, analytical);
% xlabel('X [unit]'); ylabel('P');
% legend('Numerical', 'Analytical');

