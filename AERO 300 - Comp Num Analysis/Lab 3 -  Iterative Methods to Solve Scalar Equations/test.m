%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 3 - Advance Data Loading and Plotting: 4/19/24

%% Clear Workspace

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1:

fx = @(x) x^3+1.0142*x^2-19.3629*x+15.8398;
gx = @(x) 3*x^2+2.0284*x-19.3629;
hx = @(x) 6*x+2.0824; 

d = 3; %Degree of polynomial
s = 1; %Step size
g = -6; %inital guess
tol = 1e-6;

M = bracket(g,s,d,fx);


[M1,e] = newton(d,M,fx,gx,tol);

function [roots_newton, error_newton] = newton(root_guess, roots, f, g, tol) % need to call g since newton uses it
roots_newton = zeros(root_guess,1); % stores the roots found

xO = 0;

for n = 1:root_guess
    a = roots(n,1); % defines where to place the roots
    b = roots(n,2);
    iterations = zeros(root_guess,1);

    xo = (a+b)/2;
    xoo = 0;
    while abs(f(xoo))>tol % the newton formula
        xoo = xO; % whatever xO is, set equal to xoo for last iteration
        xO = xo - f(xo)/g(xo);
        xo = xO; % redefine xo and xO to run the loop once more
        iterations(n,1) = iterations(n,1) + 1;
        w = iterations(n,1);
        iteration_newton(n,w) = iterations(n,1);
        error_newton(n,w) = xO;
    end

    roots_newton(n) = xO; % xo is redefined to be the roots 
        
end
end

%% g

function [M] = bracket(g, s, d, fx) %g is initial guess, s is step, d is # roots, fx is function
    format long
    M = zeros(d,2);  
    lBracket = g;
    d1 = 1;
    
    while d1 ~= d+1
        hBracket = lBracket+s;
        h = fx(hBracket);
        l = fx(lBracket);
        done = 0;
        if h*l<0 %if it is <0 then ans is negative and there is root in bracket
            while done == 0 
                hBracket2 = lBracket + 1e-6;
                h = fx(hBracket2);
                l = fx(lBracket);
                if h*l<0
                    M(d1,1) = lBracket;
                    M(d1,2) = hBracket2;
                    lBracket = hBracket2;
                    done = 1;
                    d1 = d1 + 1;
                else
                    lBracket = hBracket2;
                end
            end
        elseif h*l == 0
            if h == 0
                M(d,1) = h;
                M(d,2) = h;
            elseif l == 0
                M(d,1) = l;
                M(d,2) = l;
            else
                disp('error 0');
            end
        else
            lBracket = hBracket;
        end
        
    end

end






