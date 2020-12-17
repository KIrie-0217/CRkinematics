global offset
global rc
global N
global L
global h
%æ—¢å®šå?¤ã®è¨­å®?
offset = [20,60,100,0,40,80,90,130,170];
offset = offset * pi /180;
rc = [4.6,4.6,4.6,4.8,4.8,4.8,4.6,4.6,4.6];
N = 5;
L = [120,120,120,100,100,100,75,75,75];
h = [1, 1, 1,1,1,1,1,1,1];

%theta = table2array(follow2d_jun(3));
phi = [0 0 0];

theta = [45 0 0];
wire_length = inv_LSK(theta,phi)


    wire_11 =[];
    wire_12 =[];
    wire_13 =[];
    wire_21 =[];
    wire_22 =[];
    wire_23 =[];
    wire_31 =[];
    wire_32 =[];
    wire_33 =[];
    
    tsin  = [];
for i = 1:201
    if i >100
        wire_11 =[wire_11;wire_length(1,1)];
        wire_12 =[wire_12;wire_length(2,1)];
        wire_13 =[wire_13;wire_length(3,1)];
        wire_21 =[wire_21;wire_length(1,2)];
        wire_22 =[wire_22;wire_length(2,2)];
        wire_23 =[wire_23;wire_length(3,2)];
        wire_31 =[wire_31;wire_length(1,3)];
        wire_32 =[wire_32;wire_length(2,3)];
        wire_33 =[wire_33;wire_length(3,3)];
    elseif i >50
        wire_11 =[wire_11;0];
        wire_12 =[wire_12;0];
        wire_13 =[wire_13;0];
        wire_21 =[wire_21;wire_length(1,2)];
        wire_22 =[wire_22;wire_length(2,2)];
        wire_23 =[wire_23;wire_length(3,2)];
        wire_31 =[wire_31;wire_length(1,3)];
        wire_32 =[wire_32;wire_length(2,3)];
        wire_33 =[wire_33;wire_length(3,3)];
    else
        wire_11 =[wire_11;0];
        wire_12 =[wire_12;0];
        wire_13 =[wire_13;0];
        wire_21 =[wire_21;0];
        wire_22 =[wire_22;0];
        wire_23 =[wire_23;0];
        wire_31 =[wire_31;wire_length(1,3)];
        wire_32 =[wire_32;wire_length(2,3)];
        wire_33 =[wire_33;wire_length(3,3)];
    end
    tsin = [tsin;(i-1)*0.01];
end

w11 = timeseries(wire_11,tsin);
w12 = timeseries(wire_12,tsin);
w13 = timeseries(wire_13,tsin);
w21 = timeseries(wire_21,tsin);
w22 = timeseries(wire_22,tsin);
w23 = timeseries(wire_23,tsin);
w31 = timeseries(wire_31,tsin);
w32 = timeseries(wire_32,tsin);
w33 = timeseries(wire_33,tsin);


beta = fw_LSK(wire_length);
pose = pose_CR(beta);

x = pose(1,1:16);
y = pose(2,1:16);
z = pose(3,1:16);
plot3(x,y,z,'-o');
xlabel("x")
ylabel("y")
zlabel("z")
grid on 
daspect([1 1 1])
hold on 

z2 = [1:1:350];     
x2 = 50*cos(z2*2*pi()/350)-50;
y2 = z2*0;

plot3(x2,y2,z2)

for i = 1:size(wire_length,2)
    
end


hold off
