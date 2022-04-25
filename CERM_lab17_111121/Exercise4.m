% CERM_Lab_111121
% Exercise 4

clc
clear all
close all
format short

            % Defining A

A = (10 ^ -8) * eye(20) + hilb(20);
Xex = ones(20, 1);
gamma = 10 ^ -3;

e = gamma * randn(20, 1);

Xper = Xex + e;
???