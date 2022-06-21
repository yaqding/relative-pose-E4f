clear;
close all;
clc
sigma = 0; % noise level
parallax = 0.5;
fov = 90; npt = 200; 
iter = 100;

for kk = 1:iter   %
    f1 = randi([300 3000]); f2 = randi([300 3000]);
    K1 = [f1 0 0; 0 f1 0; 0 0 1]; K2 = [f2 0 0; 0 f2 0; 0 0 1]; 
    [X1, X2, tr, X0, Tx1, Tz1, Tx2, Tz2] = generate_pt(sigma, parallax,fov,npt,K1,K2);
    
    Rg = tr(1:3,1:3); % truth
    Tg = tr(1:3,4);
    
    xc = X1;
    xd = X2;
    xc(3,:) = 1; xd(3,:) = 1;
    fn = max(xd(:));
    Kn = [fn 0 0; 0 fn 0; 0 0 1];
    
    x1 = K1\xc; x2 = Kn\xd;
    Ra = (Tz1*Tx1).';
    Rb = (Tz2*Tx2).';
    
    Rg = Rb.'*Rg*Ra;
    Tg = Rb.'*Tg;
    
    
    ind = randsample(size(x1,2),4);
    
    %%  
    res1 = E4f_test(x1(1:3,ind), x2(1:3,ind), Ra, Rb, Rg, Tg,fn, f2);
    
    
end

