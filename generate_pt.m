function [im1, im2, T, x0, Tx1, Tz1, Tx2, Tz2] = generate_pt(sigma, parallax,fov,npt,K1,K2)


a = rand(1,npt)*6-3;
b = rand(1,npt)*6-3;
c = rand(1,npt)*4+2;
x0 = [a;b;c];

T = zeros(3,4);

ay = rand*20-10;
ax1 = rand*20-10; az1 = rand*20-10; ax2 = rand*20-10; az2 = rand*20-10;

T(1:3,1:3) = roty(ay);
Tx1 = rotx(ax1);
Tz1 = rotz(az1);
Tx2 = rotx(ax2);
Tz2 = rotz(az2);
A = rand(3,1);
T(1:3,4) = A/norm(A)*parallax; %



fo1 = K1(1,1)*tand(fov/2);
fo2 = K2(1,1)*tand(fov/2);


x1a = Tz1*Tx1*x0;
x2a = Tz2*Tx2*(T(:,1:3)*x0+T(:,4));

x1 = x1a./repmat(x1a(3,:),3,1);
x2 = x2a./repmat(x2a(3,:),3,1);

x1_image = K1*x1; % image points
x2_image = K2*x2;

x1_image_noi(1:2,:) = x1_image(1:2,:) + normrnd(0,sigma,2,npt);%
x2_image_noi(1:2,:) = x2_image(1:2,:) + normrnd(0,sigma,2,npt);%

% feild of view -
I = find(x1_image_noi(1,:)<fo1 & x1_image_noi(2,:)<fo1 & x1_image_noi(1,:)>-fo1 & x1_image_noi(2,:)>-fo1 ...
    & x2_image_noi(1,:)<fo2 & x2_image_noi(2,:)<fo2 & x2_image_noi(1,:)>-fo2 & x2_image_noi(2,:)>-fo2);

im1 = x1_image_noi(:,I);
im2 = x2_image_noi(:,I);




end