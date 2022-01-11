%% Averaging function. Required for mapping u values on v cell locations and vice versa for the convection operator.
function B = Averaging(A,k,nx,ny)
%k denotes whether we are averaging u or v. Input must be a row vector.
% Statu
if k == 1
 Ares = reshape(A,ny+2,[]);
 Aresn = [zeros(size(Ares,1),1),Ares,zeros(size(Ares,1),1)]; % Added additional row of zeros for averaging cells with no 'compliment'.
 AAvg1 = 0.5*(Aresn(:,1:end-1)+Aresn(:,2:end));
 AAVg2 = 0.5*(AAvg1(1:end-1,:)+AAvg1(2:end,:));
 B = reshape(AAVg2,[],1);
 sizein = size(Ares);
 sizeout = size(AAVg2);
elseif k == 2
 Ares = reshape(A,ny+1,[]);
 Aresn = [zeros(1,size(Ares,2));Ares;zeros(1,size(Ares,2))]; % Added additional column of zeros for averaging cells with no 'compliment'.
 AAvg1 = 0.5*(Aresn(1:end-1,:)+Aresn(2:end,:));
 AAVg2 = 0.5*(AAvg1(:,1:end-1)+AAvg1(:,2:end));
 B = reshape(AAVg2,[],1);
 sizein = size(Ares);
 sizeout = size(AAVg2);
end
