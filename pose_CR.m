function pose = pose_CR(beta)
    global L 
    global h
    global N
    xyz = [0;0;0;1];
    pose = [xyz];
    beta = beta/N*pi/180;
    beta = round(beta,3);
    beta = beta + 1e-5;
  
    L_ = L/10;
    simtransi = [];
    for i = 1:size(beta,2)
        
        simtrans =  Rx(beta(2,i)) * tz(L_(i))  * Ry(beta(1,i)) * tz(L_(i));
        if i == 1
            simtransi = simtrans^5;
            pose = [pose simtrans*xyz simtrans^2*xyz simtrans^3*xyz simtrans^4*xyz simtrans^5*xyz];
        else
            pose = [pose simtransi*simtrans*xyz simtransi*simtrans^2*xyz simtransi*simtrans^3*xyz ...
                simtransi*simtrans^4*xyz simtransi*simtrans^5*xyz];
            simtransi = simtransi*simtrans^5;
        %simtrans_ = ty(L(i)/beta(2,i)) * Rx(beta(2,i)) * ty(-L(i)/beta(2,i)) * tz(h(i)) ...
        %* tx(L(i)/beta(1,i)) * Ry(beta(1,i)) * tx(-L(i)/beta(1,i)) * tz(h(i)) * Rz(pi/2);
        end
        
    end
end

function a = Rx(c)
    a = [ 1 0 0 0;
          0 cos(c) -sin(c) 0;
          0 sin(c) cos(c) 0;
          0 0 0 1];
end
function a = Ry(c)
    a = [ cos(c) 0 sin(c) 0;
          0 1 0 0;
          -sin(c) 0 cos(c) 0;
          0 0 0 1];
end
function a = Rz(c)
    a = [ cos(c) -sin(c) 0 0;
          sin(c) cos(c) 0 0;
          0 0 1 0;
          0 0 0 1];
end
function a = tx(h)
    a = [ 1 0 0 h;
          0 1 0 0;
          0 0 1 0;
          0 0 0 1];
end
function a = ty(h)
    a = [ 1 0 0 0;
          0 1 0 h;
          0 0 1 0;
          0 0 0 1];
end
function a = tz(h)
    a = [1 0 0 0;
         0 1 0 0;
         0 0 1 h;
         0 0 0 1];
end


