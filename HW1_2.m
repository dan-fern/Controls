%%% Homework 1.2 uses matrix multiplication to change the basis of the
%%% provided matrices.

clear all
close all
home

beta = [-4,-2,-3;15,7,10;-5,-2,-3];

Aold = [18,9,13;50,23,35;-65,-31,-46];
Bold = [-1;0;1];
Cold = [5,-5,5];

Anew = beta * Aold * (beta .^ -1)
Bnew = beta * Bold
Cnew = Cold * beta .^ -1
