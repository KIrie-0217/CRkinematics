function pose = pose_CR2d(theta)
    global N
    xy = [0;0;1];
    pose = [xy];
    theta = table2array(theta);
    beta = theta/5*pi/180;
    beta = round(beta,3);
    beta = beta + 1e-5;


    for i = 1:numel(theta)

        simtrans =  t2d(12)  * R2d(beta(i)) * t2d(12);
        if i == 1
            simtransi = simtrans^5;
            pose = [pose t2d(12)*simtrans^2*xy  simtrans^5*xy];
        else
            pose = [pose t2d(12)*simtransi*simtrans^2*xy ...
                 simtransi*simtrans^5*xy];
            simtransi = simtransi*simtrans^5;
        %simtrans_ = ty(L(i)/beta(2,i)) * Rx(beta(2,i)) * ty(-L(i)/beta(2,i)) * tz(h(i)) ...
        %* tx(L(i)/beta(1,i)) * Ry(beta(1,i)) * tx(-L(i)/beta(1,i)) * tz(h(i)) * Rz(pi/2);
        end
    end

end