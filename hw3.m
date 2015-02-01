%%% Daniel Fernández
%%% ME531
%%% PS #3

%%  Block 3.1
clear all
close all
home

k1 = 2; k2 = 3;  % spring constants
c1 = 4; c2 = 4;  % damping constants
lambda = [-7, -7];  % eigenvalues
xin = [1/2;5/6]; 
xo = [0;0];

damping = [(c1 + c2), -c2; -c2, c2];
stiffness = [-(k1 + k2), k2; k2, -k2];
input = [0;1];
Ci = inv(damping);

A = damping \ stiffness;       % get A matrix
B = damping \ input;   % get B matrix

P = [B, A*B];
Pi = inv(P);     % i for inverse
p = Pi(2,:);
Mi = [p; p * A]; % i for inverse
M = inv(Mi);
Ac = Mi * A * M; Ac(1) = 0;  % Ac for CCF; set index 1 to 0 to clean it up
Bc = Mi * B; Bc(1) = 0;      % Bc for CCF; set index 1 to 0 to clean it up

syms Kc1 Kc2
% Use the solver function to populate Kc for CCF
var = Ac + Bc * [Kc1 Kc2]; 
S = solve(-var(2,1) == lambda .^ 2, var(2,2) == 2 * lambda);
Kc = [eval(S.Kc1) eval(S.Kc2)];
K = Kc * Mi;      % get gains!!

Aeq = A + B * K;  % get Aeq!!

x = []; 
% Use this loop to populate all values of x for the x_in
for i = 0.01:0.01:2.5
    xi = expm(Aeq * i) * xo - Aeq ^ (-1) * expm(Aeq * i) * B * K * xin ...
        + Aeq ^ (-1) * B * K * xin;
    x = [x xi]; 
end
figure(1)
t = 0.01:0.01:2.5;
% Lastly, the plots.
plot(t, x(1,:), ':r');
hold on
plot(t, x(2,:), '-b');
xlabel('Time, s'); ylabel('Amplitude'); title('Problem 3.1');
legend('Z_1', 'Z_2', 'Location', 'NorthEast');

%%  Block 3.2
clear all
close all
home

k1 = 2; k2 = 3;  % spring constants
c1 = 4; c2 = 4;  % damping constants
lambda = [-7, -7, -7];  % eigenvalues
xin = [0;1/2;5/6]; 
xo = [0;0;0];

damping = [(c1 + c2), -c2; -c2, c2];
stiffness = [-(k1 + k2), k2; k2, -k2];
input = [0;1];
Ci = inv(damping);

Atemp = damping \ stiffness;       
Btemp = damping \ input;   

A = [0, 1, 1; 0, Atemp(1), Atemp(3); 0, Atemp(2), Atemp(4)]; % get A matrix
B = [0; Btemp(1); Btemp(2)];  % get B matrix

P = [B, A*B, A*A*B];
Pi = inv(P);     % i for inverse
p = Pi(3,:);
Mi = [p; p*A ; p*A*A]; % i for inverse
M = inv(Mi); M(7) = 0;        % set index 7 to 0 to clean it up
Ac = Mi * A * M; %Ac(1) = 0;  % Ac for CCF; set index 1 to 0 to clean it up
Bc = Mi * B; Bc(2) = 0;       % Bc for CCF; set index 2 to 0 to clean it up

syms Kc1 Kc2 Kc3
% Use the solver function to populate Kc for CCF
var = Ac + Bc * [Kc1 Kc2 Kc3]; 
S = solve(var(3,1) == lambda .^ 3, var(3,2) == -3 * lambda .^ 2, var(3,3) == 3 * lambda);
Kc = [eval(S.Kc1) eval(S.Kc2) eval(S.Kc3)];
K = Kc * Mi;      % get gains!!

Aeq = A + B * K;  % get Aeq!!

x = []; 
E = [0 -1 -1; 0 0 0; 0 0 0];  % error matrix!!
% Use this loop to populate all values of x for the x_in
for i = 0.01:0.01:2.5
    xi = expm(Aeq * i) * xo - Aeq ^ (-1) * expm(Aeq * i) * (B * K - E) * xin ...
        + Aeq ^ (-1) * (B * K - E) * xin;
    x = [x xi]; 
end
figure(2)
t = 0.01:0.01:2.5;
% Lastly, the plots.
plot(t, x(1,:), ':m');
hold on
plot(t, x(2,:), '-b');
hold on
plot(t, x(3,:), '-.r');
xlabel('Time, s'); ylabel('Amplitude'); title('Problem 3.2');
legend('Error', 'Z_1', 'Z_2', 'Location', 'SouthEast');

%%  Block 3.3
clearvars -except K
close all
home

K = [-1301, 625];
k1 = 2; k2 = 3;  % spring constants
c1 = 4; c2 = 4;  % damping constants
lambda = [-70 + 70i, -70 - 70i];  % eigenvalues
xin = [1/2;5/6]; 
xo = [0;0];

damping = [(c1 + c2), -c2; -c2, c2];
stiffness = [-(k1 + k2), k2; k2, -k2];
input = [0;1];
dampinv = inv(damping);

A = damping \ stiffness;       % get A matrix
B = damping \ input;           % get B matrix
C = [0 1];  

Q = [C; C*A];
Qi = inv(Q);     % i for inverse
q = Qi(:,2);
U = [A*q, q]; % i for inverse
Ui = inv(U);
Ao = Ui * A * U; %Ac(1) = 0;  % Ac for CCF; set index 1 to 0 to clean it up
Bo = Ui * B; %Bc(1) = 0;      % Bc for CCF; set index 1 to 0 to clean it up
Co = C * U;

syms Lo1 Lo2
% Use the solver function to populate Kc for CCF
pnom = poly([-70+70i, -70-70i]);
var = Ao + [Lo1;Lo2] * Co;
S = solve(-var(1) == pnom(2), -var(2) == pnom(3));
Lo = [eval(S.Lo1); eval(S.Lo2)];
L = U * Lo;      % get gains!!

Aeq = A + L * C;  % get Aeq!!

Xdot = @(T,X) A * X + B * K * (xin - X) - L * C * (xin - X);

[t, x] = ode45(Xdot, [0 .15], xo);

figure(3)
plot(t, x(:,1),'b')
hold on
plot(t, x(:,2),'r')
legend('Z_1', 'Z_2'); axis([0 0.1 -6 8]);
xlabel('Time, s'); ylabel('Amplitude'); title('Problem 3.3c');

xo = [0;0;2;0]; % now with the observer

Xdot2 = @(T,X) [A * X(1:2) + B * -K * (xin - X(3:4)); ...
    A * X(3:4) - B * K * (xin - X(3:4)) - L * (C * X(1:2) - C * X(3:4))];

[t, x] = ode45(Xdot2, [0 1], xo);

figure(4)
count = 1;
plot(t(1:count:end), x(1:count:end,1),'b'); hold on;
plot(t(1:count:end), x(1:count:end,2),'r');
plot(t(1:count:end), x(1:count:end,3),'m');
plot(t(1:count:end), x(1:count:end,4),'c');
legend('Z_1', 'Z_2', 'Observed_1', 'Observed_2', 'Location', 'SouthEast');
xlabel('Time, s'); ylabel('Amplitude'); title('Problem 3.3d');

figure(5)
count = 1;
plot(t(1:count:end), x(1:count:end,1),'b'); hold on;
plot(t(1:count:end), x(1:count:end,2),'r');
plot(t(1:count:end), x(1:count:end,3),'m');
plot(t(1:count:end), x(1:count:end,4),'c');
legend('Z_1', 'Z_2', 'Observed_1', 'Observed_2', 'Location', 'SouthEast');
axis([0 0.25 -15 10]);
xlabel('Time, s'); ylabel('Amplitude'); title('Problem 3.3d, Enhanced');



























