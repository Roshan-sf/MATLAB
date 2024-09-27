%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 4 - Iterative Methods to Solve Matrix Equations: 4/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 



% 
% syms y
% 
% solve('x+y=2',y)
% 
% 
% A = [3 -1 2 ; 1 4 2 ; -1 3 5];
% 
% 
% 
% u = s(1,1);
% ig = zeros(u,1);
% 
% [L, U] = lu(A);

% Gauss-Seidel method

 % n=input('Enter number of equations, n:  ');
 % A = zeros(n,n+1);
 % x1 = zeros(n);
 % tol = input('Enter the tolerance, tol: ');
 % m = input('Enter maximum number of iterations, m:  ');
 % 
 % A=[4 2 3 8; 3 -5 2 -14; -2 3 8 27];
 % x1=[0 0 0];
 % 
 % k = 1;
 % while  k <= m
 %    err = 0;
 %    for i = 1 : n 
 %       s = 0;
 %       for j = 1 : n 
 %          s = s-A(i,j)*x1(j);
 %       end
 %       s = (s+A(i,n+1))/A(i,i);
 %       if abs(s) > err 
 %         err  = abs(s);
 %       end
 %       x1(i) = x1(i) + s;
 %    end
 % 
 %    if err <= tol 
 %      break;
 %    else
 %      k = k+1;
 %    end
 % end 
 % 
 % fprintf('Solution vector after %d iterations is :\n', k-1);
 % for i = 1 : n 
 %   fprintf(' %11.8f \n', x1(i));
 % end






A = [3 -1 2 ; 1 4 2 ; -1 3 5];
b = [1;1;1];

z = cat(2,A,b);



 x=[0 0 0];
 
 k = 1;
 while  k <= m
    err = 0;
    for i = 1 : n 
       x1 = 0;
       for j = 1 : n 
          x1 = x1-A2(i,j)*x(j);
       end
       x1 = (x1+A2(i,n+1))/A2(i,i);
       if abs(x1) > err 
         err  = abs(x1);
       end
       x(i) = x(i) + x1;
    end

    if err <= tol 
      break;
    else
      k = k+1;
    end
 end 























 n=s(1,1);
 m = 10000;
 A2 = cat(2,A,b); %Combining to single matrix for easier calcs

 x1=[0 0 0];
 
 k = 1;
 while  k <= m
    err = 0;
    for i = 1 : n 
       s = 0;
       for j = 1 : n 
          s = s-A2(i,j)*x1(j);
       end
       s = (s+A2(i,n+1))/A2(i,i);
       if abs(s) > err 
         err  = abs(s);
       end
       x1(i) = x1(i) + s;
    end

    if err <= tol 
      break;
    else
      k = k+1;
    end
 end 
 
 fprintf('Solution vector after %d iterations is :\n', k-1);
 for i = 1 : n 
   fprintf(' %11.8f \n', x1(i));
 end









