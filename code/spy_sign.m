function []=spy_sign(A)
% author: Marc Hesse
% date: 12 July 2013
% This function plots the sparsity pattern of a function and indicates the
% sign by color. Green is positive and red is negative
spy(A+abs(A),'go'), hold on
spy(A-abs(A),'ro')
legend('+','-')