%% Averaging function. Required for mapping phi values on u and v cell locations.
function B = Averaging_phi(A,k,nx,ny)
if k == 1
 Ares = reshape(A,[],nx);
 %Aresn = [zeros(1,size(Ares,2));Ares;zeros(1,size(Ares,2))]; % Added additional row of zeros for averaging cells with no 'compliment'.
 AAvg1 = 0.5*(Ares(:,1:end-1)+Ares(:,2:end));
 sizein = size(Ares);
 sizeout = size(AAvg1);
 B = AAvg1;
elseif k == 2
 Ares = reshape(A,ny,[]);
 %Aresn = [zeros(size(Ares,1),1),Ares,zeros(size(Ares,1),1)]; % Added additional column of zeros for averaging cells with no 'compliment'.
 AAvg1 = 0.5*(Ares(1:end-1,:)+Ares(2:end,:));
 sizein = size(Ares);
 sizeout = size(AAvg1);
 B = AAvg1;
end