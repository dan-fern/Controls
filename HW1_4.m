clear all
close all
home

% set up the grid space
[x1, x2] = meshgrid(-5:0.5:5,-5:.5:5); 

% initialize the matrices
A{1} = [-8,-6;0,-2];
A{2} = [-8,-6;0,2];
A{3} = [0,4;-4,0];
A{4} = [4,-4;4,4];
A{5} = [1,0;0,0];
A{6} = [1,0;1,0];
A{7} = [0,0;0,0];
A{8} = [-1,1,0;0,-1,0;0,0,-1];

% creates 2 xdot functions and charts through the quiver() function.
% could not successfully troubleshoot the ode45 function, commented it out
count = 1;
for i = 1:length(A)-1
    Anew = A{i};
    x1dot = Anew(1) * x1 + Anew(3) * x2;
    x2dot = Anew(2) * x1 + Anew(4) * x2;
    figure(count);
    quiver(x1, x2, x1dot, x2dot);
%     hold on
%     IC = [1,1;1,2;1,3;1,4];
%     for j = 1:length(IC(:,1))
%         [~,X] = ode45(@my_fun,[0 50],IC(j,:));
%         u = X(:,1);
%         w = X(:,2);
%         plot(u,w,'r')
%     end
    xlabel('x1'); ylabel('x2');
    xlim([x1(1),x1(end)]); ylim([x1(1),x1(end)])
    count = count + 1;
end

% arrows are indicative of vector direction and magnitude.
% lines are indicative of overall trajectory
% see below for 6.2.h
[x1, x2, x3] = meshgrid(-5:0.5:5,-5:.5:5,-5:.5:5); 
Anew = A{8};
x1dot = Anew(1) * x1 + Anew(4) * x2 + Anew(7) * x3;
x2dot = Anew(2) * x1 + Anew(5) * x2 + Anew(8) * x3;
x3dot = Anew(3) * x1 + Anew(6) * x2 + Anew(9) * x3;
figure(8);
quiver3(x1,x2,x3,x1dot,x2dot,x3dot);
xlabel('x1'); ylabel('x2'); zlabel('x3');
xlim([x1(1),x1(end)]); ylim([x1(1),x1(end)]); zlim([x1(1),x1(end)]);

clear i;
clear j;
clear count;
