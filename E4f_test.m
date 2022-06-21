function res = E4f_test(p1, p2, g1, g2, Rg, Tg, fn, f2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

q1 = g1*p1;
q1 = q1./q1(3,:);
a1 = g2(1,2)/g2(1,1);
a2 = g2(2,1)/g2(1,1);
a3 = g2(2,2)/g2(1,1);
a4 = g2(2,3)/g2(1,1);
a5 = g2(3,1)/g2(1,1);
a6 = g2(3,2)/g2(1,1);
a7 = g2(3,3)/g2(1,1);
q2 = [1 a1; a2 a3; a5 a6]*p2(1:2,:);
data = [a1 a2 a3 a4 a5 a6 a7 q1(1,:) q1(2,:) q2(1,:) q2(2,:) q2(3,:)];
sols = solver_e4f_polyeig(data);

if (size(sols,2) >0)
    
    for k = 1:size(sols,2)
        R(:,:,k) = [1-sols(1,k)^2 0 2*sols(1,k); 0 1+sols(1,k)^2 0; -2*sols(1,k) 0 1-sols(1,k)^2];
        
        Re = (g2.')*R(:,:,k)*g1/(1+sols(1,k)^2);
        
        A(1,:,k) = cross([p2(1,1) p2(2,1) sols(2,k)]*(g2.'),(R(:,:,k)*g1*[p1(1,1) p1(2,1) 1].').');
        A(2,:,k) = cross([p2(1,2) p2(2,2) sols(2,k)]*(g2.'),(R(:,:,k)*g1*[p1(1,2) p1(2,2) 1].').');
        A(3,:,k) = cross([p2(1,3) p2(2,3) sols(2,k)]*(g2.'),(R(:,:,k)*g1*[p1(1,3) p1(2,3) 1].').');
        error_r(1,k) = acosd((trace(Rg*(Re.'))-1)/2);
        error_f(1,k) = abs(f2-fn*sols(2,k))/f2;
        try
            [~,~,V] = svd(A(:,:,k),0);
            t(:,k) = (g2.')*V(:,3);
            error_t(1,k) = acosd(  abs(Tg.'*t(:,k))/(norm(Tg)*norm(t(:,k))));
        catch
            continue;
        end
        
    end
    
    res = [ error_r; error_t; error_f];
    
else
    
    res = [99 99 99]';
    
end



end

