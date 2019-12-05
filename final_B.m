% This is for Part B
N = 3;
Ae = zeros(N,5);
Be = zeros(N,5);
Ao = zeros(N,5);
Bo = zeros(N,5);

for n = 0:1:N-1
    P = 3;
%     n = 0;
    syms a q a1 a2 a3 a4 a5;
    a = -(n^2 + a1*q + a2*q^2 + a3*q^3 + a4*q^4 + a5*q^5);
    Aeven = sym(zeros(P));
    Beven = sym(zeros(P));
    Aodd = sym(zeros(P));
    Bodd = sym(zeros(P));
    Aeven(1, 1:2) = [-a, q];
    Beven(1, 1:2) = [-a-4, q];
    Aodd(1, 1:2) = [-a-1+q, q];
    Bodd(1, 1:2) = [-a-1-q, q];
    for i = 2:1:P-1
        if i == 2
            Aeven(i, i-1:i+1) = [2*q, -a-(2*i-2)^2, q];
        else
            Aeven(i, i-1:i+1) = [q, -a-(2*i-2)^2, q];
        end
        Beven(i, i-1:i+1) = [q, -a-(2*i)^2, q];
        Aodd(i, i-1:i+1) = [q, -a-(2*i-1)^2, q];
        Bodd(i, i-1:i+1) = [q, -a-(2*i-1)^2, q];
    end
    Aeven(P, P-1:P) = [q, -a-(2*P-2)^2];
    Beven(P, P-1:P) = [q, -a-(2*P)^2];
    Aodd(P, P-1:P) = [q, -a-(2*P-1)^2];
    Bodd(P, P-1:P) = [q, -a-(2*P-1)^2];

    detAeven = det(Aeven);
    disp(detAeven);
    detBeven = det(Beven);
    disp(detBeven);
    detAodd = det(Aodd);
    disp(detAodd);
    detBodd = det(Bodd);
    disp(detBodd);
    
    [c,~]= poly_coeffs(detAeven,'q');
    [ae10,ae20,ae30,ae40,ae50] = solve(c(end-5),c(end-4),c(end-3),c(end-2),c(end-1),a1,a2,a3,a4,a5);
    Ae(n+1, 1) = ae10;
    Ae(n+1, 2) = ae20;
    Ae(n+1, 3) = ae30;
    Ae(n+1, 4) = ae40;
    Ae(n+1, 5) = ae50;
    
    [c,~]= poly_coeffs(detBeven,'q');
    [be10,be20,be30,be40,be50] = solve(c(end-5),c(end-4),c(end-3),c(end-2),c(end-1),a1,a2,a3,a4,a5);
    Be(n+1, 1) = be10;
    Be(n+1, 2) = be20;
    Be(n+1, 3) = be30;
    Be(n+1, 4) = be40;
    Be(n+1, 5) = be50;
    
    [c,~]= poly_coeffs(detAodd,'q');
    [ao10,ao20,ao30,ao40,ao50] = solve(c(end-5),c(end-4),c(end-3),c(end-2),c(end-1),a1,a2,a3,a4,a5);
    Ao(n+1, 1) = ao10;
    Ao(n+1, 2) = ao20;
    Ao(n+1, 3) = ao30;
    Ao(n+1, 4) = ao40;
    Ao(n+1, 5) = ao50;
    
    [c,~]= poly_coeffs(detBodd,'q');
    [bo10,bo20,bo30,bo40,bo50] = solve(c(end-5),c(end-4),c(end-3),c(end-2),c(end-1),a1,a2,a3,a4,a5);
    Bo(n+1, 1) = bo10;
    Bo(n+1, 2) = bo20;simp
    Bo(n+1, 3) = bo30;
    Bo(n+1, 4) = bo40;
    Bo(n+1, 5) = bo50;
end



x = -5:0.05:5;
y1 = -(0 + x * Ae(1,1) + Ae(1,2) * x.*x + Ae(1,3) * x.*x.*x + Ae(1,4) * x.*x.*x.*x + Ae(1,5) * x.*x.*x.*x.*x);
% y2 = -(1 + x * Ae(2,1) + Ae(2,2) * x.*x + Ae(2,3) * x.*x.*x + Ae(2,4) * x.*x.*x.*x + Ae(2,5) * x.*x.*x.*x.*x);
y3 = -(4 + x * Ae(3,1) + Ae(3,2) * x.*x + Ae(3,3) * x.*x.*x + Ae(3,4) * x.*x.*x.*x + Ae(3,5) * x.*x.*x.*x.*x);

% y4 = -(0 + x * Be(1,1) + Be(1,2) * x.*x + Be(1,3) * x.*x.*x + Be(1,4) * x.*x.*x.*x + Be(1,5) * x.*x.*x.*x.*x);
% y5 = -(1 + x * Be(2,1) + Be(2,2) * x.*x + Be(2,3) * x.*x.*x + Be(2,4) * x.*x.*x.*x + Be(2,5) * x.*x.*x.*x.*x);
y6 = -(4 + x * Be(3,1) + Be(3,2) * x.*x + Be(3,3) * x.*x.*x + Be(3,4) * x.*x.*x.*x + Be(3,5) * x.*x.*x.*x.*x);

% y7 = -(0 + x * Ao(1,1) + Ao(1,2) * x.*x + Ao(1,3) * x.*x.*x + Ao(1,4) * x.*x.*x.*x + Ao(1,5) * x.*x.*x.*x.*x);
y8 = -(1 + x * Ao(2,1) + Ao(2,2) * x.*x + Ao(2,3) * x.*x.*x + Ao(2,4) * x.*x.*x.*x + Ao(2,5) * x.*x.*x.*x.*x);
% y9 = -(4 + x * Ao(3,1) + Ao(3,2) * x.*x + Ao(3,3) * x.*x.*x + Ao(3,4) * x.*x.*x.*x + Ao(3,5) * x.*x.*x.*x.*x);

% y10 = -(0 + x * Bo(1,1) + Bo(1,2) * x.*x + Bo(1,3) * x.*x.*x + Bo(1,4) * x.*x.*x.*x + Bo(1,5) * x.*x.*x.*x.*x);
y11 = -(1 + x * Bo(2,1) + Bo(2,2) * x.*x + Bo(2,3) * x.*x.*x + Bo(2,4) * x.*x.*x.*x + Bo(2,5) * x.*x.*x.*x.*x);
% y12 = -(4 + x * Bo(3,1) + Bo(3,2) * x.*x + Bo(3,3) * x.*x.*x + Bo(3,4) * x.*x.*x.*x + Bo(3,5) * x.*x.*x.*x.*x);

plot(y1, x, y3, x, y6, x, y8, x, y11, x);
title("Transition curves in Mathieu's equation");
axis([-6, 2, -5, 5]);
xlabel('a');
ylabel('q');
% P = 3;
% n = 1;
% syms a q a1 a2 a3 a4 a5;
% %a = n^2/4 + a1*q + a2*q^2 + a3*q^3 + a4*q^4 + a5*q^5;
% Aodd = sym(zeros(P));
% Aodd(1, 1:2) = [a-1/4+q/2, q/2];
% for i = 2:1:P-1
%     Aodd(i, i-1:i+1) = [q/2, a-(2*i-1)^2/4, q/2];
% end
% Aodd(P, P-1:P) = [q/2, a-(2*P-1)^2/4];
% Aodd
% detA = det(Aodd);
% disp(detA);
% %fa = fcontour(detA);
% [c,t]= poly_coeffs(detA,'q');
% [a10,a20,a30,a40,a50] = solve(c(end-5),c(end-4),c(end-3),c(end-2),c(end-1),a1,a2,a3,a4,a5);
% a10
% a20
% a30
% a40
% a50