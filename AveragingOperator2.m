%% Averaging function. Required for mapping u values on cell centre locations for the gradient operator.
%Status - COMPLETE  
function [B ,sizein, sizeout] = AveragingOperator2(A,k,nx,ny)
if k == 1
 Ares = reshape(A,[],nx+1);
 Aresn = [zeros(1,size(Ares,2));Ares;zeros(1,size(Ares,2))]; % Added additional row of zeros for averaging cells with no 'compliment'.
 AAvg1 = 0.5*(Aresn(1:end-1,:)+Aresn(2:end,:));
 B = reshape(AAvg1,[],1);
 sizein = size(Ares);
 sizeout = size(AAvg1);
elseif k == 2
 Ares = reshape(A,ny+1,[]);
 Aresn = [zeros(size(Ares,1),1),Ares,zeros(size(Ares,1),1)]; % Added additional column of zeros for averaging cells with no 'compliment'.
 AAvg1 = 0.5*(Aresn(:,1:end-1)+Aresn(:,2:end));
 B = reshape(AAvg1,[],1);
 sizein = size(Ares);
 sizeout = size(AAvg1);
end