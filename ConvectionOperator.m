%% Function to compute the convection operator H. The k denotes whether the operator is being run on the x NS equation or the y NS equation.
% Status- COMPLETE
function [H] = ConvectionOperator(u,v,k,dx,dy,nx,ny)
% c denotes the component of velocity tested.
% Averaging using linear interpolation (u and v)
%uavg = 
H1 = zeros(size(u,1),size(u,1));
H2 = zeros(size(u,1),size(u,1));

if k == 1
for i = 4:size(u,1)-4
    if u(i) >=0
        % uux components
        H1(i,i+2) = -6/120*dx;
        H1(i,i+1) = 60/120*dx;
        H1(i,i) = 40/120*dx;
        H1(i,i-1) = -120/120*dx;
        H1(i,i-2) = -30/120*dx;
        H1(i,i-3) = -4/120*dx;
        
        % vuy components
        H2(i+2,i) = -6/120*dy;
        H2(i+1,i) = 60/120*dy;
        H2(i,i) = 40/120*dy;
        H2(i-1,i) = -120/120*dy;
        H2(i-2,i) = -30/120*dy;
        H2(i-3,i) = -4/120*dy;
   elseif u(i) < 0
        % uux components
        H1(i+2) = -6/120*dx;
        H1(i+1) = 60/120*dx;
        H1(i,i) = 40/120*dx;
        H1(i,i-1) = -120/120*dx;
        H1(i,i-2) = -30/120*dx;
        H1(i,i-3) = -4/120*dx;

        %vuy components
        H1(i+2) = -6/120*dy;
        H1(i+1) = 60/120*dy;
        H1(i,i) = 40/120*dy;
        H1(i,i-1) = -120/120*dy;
        H1(i,i-2) = -30/120*dy;
        H1(i,i-3) = -4/120*dy;
    end  
        % Treating interior cells adjacent to boundaries H1
        for k = nx+3:nx+1:(nx+1)*(ny+1)+2
         if u(k) >= 0
            H1(k,k) = 0;
            H1(k,k-1) = -1/(2*dx);
            H1(k,k+1) =  1/(2*dx);
         elseif u(k) < 0
            H1(k,k-1) = -2/(6*dx);
            H1(k,k) = -3/(6*dx);
            H1(k,k+1) =  6/(6*dx);
            H1(k,k+2) = -1/(6*dx);
         end
        end
        for k = 2*(nx+1)-1:nx+1:(nx+1)*(ny+1)-1
         if u(k) >= 0
            H1(k,k) = 0;
            H1(k,k-1) = -1/(2*dx);
            H1(k,k+1) =  1/(2*dx);
         elseif u(k) < 0
            H1(k,k+1) = 2/(6*dx);
            H1(k,k) = -3/(6*dx);
            H1(k,k-1) = -6/(6*dx);
            H1(k,k-2) =  1/(6*dx);
         end
        end

    % Treating ghost cells H1
    H1(1:nx+1:(nx+1)*(ny+1)+1,:) = 0;
    H1((nx+1):nx+1:(nx+1)*(ny+2),:) = 0;
    H1(1:nx+1,:) = 0;
    H1((nx+1)*(ny+1)+1:(nx+1)*(ny+2),:) = 0;

    % Treating interior cellls adjacent to the boundary H2
    for k1 = nx+3:1:2*(nx+1)
        if u(k) >= 0
            H2(k,k) = 0;
            H2(k-1,k) = -1/(2*dx);
            H2(k+1,k) =  1/(2*dx);
         elseif u(2) < 0
            H2(k-1,k) = -2/(6*dx);
            H2(k,k) = -3/(6*dx);
            H2(k+1,k) =  6/(6*dx);
            H2(k+2,k) = -1/(6*dx);
        end
    end
    for k = (nx+1)*(ny)+1:1:(nx+1)*(ny+1)
         if u(k) > 0
            H2(k,k) = 0;
            H2(k-1,k) = -1/(2*dx);
            H2(k-1,k) =  1/(2*dx);
         elseif u(k) < 0
            H2(k+1,k) = 2/(6*dx);
            H2(k,k) = -3/(6*dx);
            H2(k-1,k) = -6/(6*dx);
            H2(k-2,k) =  1/(6*dx);
         end
    end
     % Treating Ghost Cells H2
     H2(1:nx+1:(nx+1)*(ny+1)+1,:) = 0;
     H2((nx+1):nx+1:(nx+1)*(ny+2),:) = 0;
     H2(1:nx+1,:) = 0;
     H2((nx+1)*(ny+1)+1:(nx+1)*(ny+2),:) = 0;

     % Deleting invalid cells
        H2([1 nx+1 ((nx+1)*(ny+1)+1) (ny+2)*(nx+1)],:) = 0;
        H1([1 nx+1 ((nx+1)*(ny+1)+1) (ny+2)*(nx+1)],:) = 0;   
end
size(H1);
size(H2);
B = Averaging(v,2,nx,ny);
 %H = 
 H = u.*(H1*u)+B.*(H2*u);
 %H = H1;

elseif k == 2
    for i = 4:size(v,1)-4
      if abs(v(i))/v(i) > 0
        %uvx components
        H1(i,i+2) = -6/120*dx;
        H1(i,i+1) = 60/120*dx;
        H1(i,i) = 40/120*dx;
        H1(i,i-1) = -120/120*dx;
        H1(i,i-2) = -30/120*dx;
        H1(i,i-3) = -4/120*dx;
        
        % vvy components
        H2(i+2,i) = -6/120*dy;
        H2(i+1,i) = 60/120*dy;
        H2(i,i) = 40/120*dy;
        H2(i-1,i) = -120/120*dy;
        H2(i-2,i) = -30/120*dy;
        H2(i-3,i) = -4/120*dy;
      elseif abs(u(i))/u(i) < 0
        % uvx components
        H1(i+2) = -6/120*dx;
        H1(i+1) = 60/120*dx;
        H1(i,i) = 40/120*dx;
        H1(i,i-1) = -120/120*dx;
        H1(i,i-2) = -30/120*dx;
        H1(i,i-3) = -4/120*dx;

        %vvy components
        H2(i+2) = -6/120*dy;
        H2(i+1) = 60/120*dy;
        H2(i,i) = 40/120*dy;
        H2(i,i-1) = -120/120*dy;
        H2(i,i-2) = -30/120*dy;
        H2(i,i-3) = -4/120*dy;
      end  
    end
       % Treating interior cells adjacent to boundaries H1
        for k = nx+4:nx+2:(nx+2)*(ny-1)+2
         if u(k) > 0
            H1(k,k) = 0;
            H1(k,k-1) = -1/(2*dx);
            H1(k,k+1) =  1/(2*dx);
         elseif u(k) < 0
            H1(k,k-1) = -2/(6*dx);
            H1(k,k) = -3/(6*dx);
            H1(k,k+1) =  6/(6*dx);
            H1(k,k+2) = -1/(6*dx);
         end
        end
        for k = 2*(nx+2)-1:nx+2:(nx+2)*(ny)-1
         if u(k) > 0
            H1(k,k) = 0;
            H1(k,k-1) = -1/(2*dx);
            H1(k,k+1) =  1/(2*dx);
         elseif u(k) < 0
            H1(k,k+1) = 2/(6*dx);
            H1(k,k) = -3/(6*dx);
            H1(k,k-1) = -6/(6*dx);
            H1(k,k-2) =  1/(6*dx);
         end
        end
         % Treating ghost cells H1
         H1(1:nx+2:(nx+2)*(ny)+1,:) = 0;
         H1((nx+2):nx+2:(nx+2)*(ny+1),:) = 0;
         H1(1:nx+2,:) = 0;
         H1((nx+2)*(ny)+1:(nx+2)*(ny+1),:) = 0;

    % Treating interior cellls adjacent to the boundary H2
    for k1 = nx+4:1:2*(nx+2)-1
        if u(k) > 0
            H2(k,k) = 0;
            H2(k-1,k) = -1/(2*dx);
            H2(k+1,k) =  1/(2*dx);
         elseif u(2) < 0
            H2(k-1,k) = -2/(6*dx);
            H2(k,k) = -3/(6*dx);
            H2(k+1,k) =  6/(6*dx);
            H2(k+2,k) = -1/(6*dx);
        end
    end
    for k = (nx+2)*(ny-1)+2:1:(nx+1)*(ny-1)+nx+1
         if u(k) > 0
            H2(k,k) = 0;
            H2(k-1,k) = -1/(2*dx);
            H2(k-1,k) =  1/(2*dx);
         elseif u(k) < 0
            H2(k+1,k) = 2/(6*dx);
            H2(k,k) = -3/(6*dx);
            H2(k-1,k) = -6/(6*dx);
            H2(k-2,k) =  1/(6*dx);
         end
    end
     % Treating Ghost Cells H2
     H2(1:nx+2:(nx+2)*(ny)+1,:) = 0;
     H2((nx+2):nx+2:(nx+2)*(ny+1),:) = 0;
     H2(1:nx+2,:) = 0;
     H2((nx+2)*(ny)+1:(nx+2)*(ny+1),:) = 0;

     % Deleting invalid cells
     H2([1 nx+2 ((nx+2)*(ny)+1) (ny+1)*(nx+2)],:) = 0;
     H1([1 nx+2 ((nx+2)*(ny)+1) (ny+1)*(nx+2)],:) = 0; 

     sz1 = size(H1);
     sz2 = size(H2);
     H = v.*(H2*v)+Averaging(u,1,nx,ny).*(H1*v);
end

end

