% Davis Sok
% Aero 300 LAB1 : 2024.04.24

close all;              % clears all
clear all;              % clears workspace
clc;                    % clears command line

% Please note: labeling this MATLAB assignment in chronological order is
% neigh immpossible due to MATLAB seething at the meer sight of a
% statement below a function. 

%% Presets & Calling Functions 

% Define functions
f = @(x) x^3 + 1.0142*x^2 - 19.3629*x + 15.8398; % function of choice
g = @(x) 3*x^2 + 2.0284*x - 19.3629; % first derivative 
h = @(x) 6*x+ 2.0824; % second derivative

% Define settings
initial_guess = 0; % where to start bracketing
root_guess = 3; % how many roots?
step_size = pi/10; % how much space do you want to check?
tol = 10^-6; % tolerance

% Call Bracketing Function
roots = bracketing(initial_guess, root_guess, step_size, f); % calling bracketing function ahead
disp("Approximate Root Locations:")
disp("     (a)       (b)")
disp(roots)
disp('If zeros have been returned, please choose a better initial guess or root guess.')

% Call Bisection Method Functon
for i = 1:10000 % iterating 10000 times since question asks for it
    tic
    [roots_bisection, iteration_bisection, error_bisection] = bisection(root_guess, roots, f, tol); % calling bisecting function
    time_bisection = toc;
end

time_average_bisection = time_bisection/10000; % average is the total time divided by the total iterations
disp('--------------------')
disp("Bisection Method Roots:")
disp(roots_bisection)
disp(['The average time for the Bisection Method is ' num2str(time_average_bisection) ' seconds.']);

absolute_bisection = abs(roots_bisection - error_bisection);

% Call Newton's Method Function
for i = 1:10000
    tic
    [roots_newton, iteration_newton, error_newton] = newton(root_guess, roots, f, g, tol);
    time_newton = toc;
end

time_average_newton = time_newton/10000;
disp('--------------------')
disp("Newton Method Roots:")
disp(roots_bisection)
disp(['The average time for the Newton Method is ' num2str(time_average_newton) ' seconds.']);

absolute_newton = abs(roots_newton - error_newton);

% Call Halley's Method Function
for i = 1:10000
    tic
    [roots_halley, iteration_halley, error_halley] = halley(root_guess, roots, f, g, h, tol);
    time_halley = toc;
end

time_average_halley = time_halley/10000;
disp('--------------------')
disp("Halley Method Roots:")
disp(roots_bisection)
disp(['The average time for the Halley Method is ' num2str(time_average_halley) ' seconds.']);

absolute_halley = abs(roots_halley - error_halley);

% Plot p(x) - x
x_domain = linspace(-10,8,1000); % setting up the x-domain
px = (x_domain).^3 + 1.0142*(x_domain).^2 - 19.3629*(x_domain) + 15.8398; % plugging in x-domain values into y function

figure
plot(x_domain, px)
grid on;
xlabel('x');
ylabel('p(x)');
legend('p(x)');
title('Figure 1 - x vs. p(x)');

% Plot all methods
figure
semilogy(iteration_bisection(1,:),absolute_bisection(1,:), '*r')
xlabel('Iteration');
ylabel('Absolute Error');
title('Convergence of Each Method');
grid on;
hold on;
semilogy(iteration_newton(1,:),absolute_newton(1,:), '*k')
semilogy(iteration_halley(1,:),absolute_halley(1,:), '*b')
legend('Bisection', 'Newton', 'Halley')

% Time Discussion:
% The time for the three methods are fast in general. However, the newton
% and halley methods are about 2x to 3x faster than the bisection method.
% In theory, the newton method should be slower than the halley method
% since the halley method converges faster, but (likely) due to the slight
% differences in code between the newton and halley methods, the newton can
% be faster. When considerating the overall picture, the halley and newton
% methods are essentially identical. 

% Slope Discussion:
% From the semilogy, we see that the newton and halley methods have a much
% steeper negative slope than the bisection method. Considering that the
% bisection method takes the longest to find the roots, this would mean
% that the slope directly coorrelates to how fast a method would converge
% onto its roots. The steeper the negative slope is, the faster a method
% converges. 

%% Bracketing Function

function roots = bracketing(initial_guess, root_guess, step_size, f) % bracketing function that uses the intial settings

% Test to the to right of initial guess
roots = zeros(root_guess, 2); % stores the roots in an array
a = initial_guess - 100*pi; % just so no-one screws up a guess
b = a + step_size; % checks a + whatever step-size user chose
n = 0; % start n at zero
nmax = 10000; % stop after 10000 checks

while n < root_guess
    if f(a)*f(b) < 0
        roots(n + 1,1) = b; % second column
        roots(n + 1,2) = a; % first column
        n = n + 1;
    end
    a = b; % set a = b
    b = b + step_size; % then have b iterative itself
    if b > (initial_guess + nmax)*step_size % root has been found
        break
    end
end

% Test to the to left of initial guess
a = initial_guess; % same thing but in reverse to the left
b = a - step_size;
n = 0;
nmax = 10000;

while n < root_guess
    if f(a)*f(b) < 0
        roots(n + 1,1) = b;
        roots(n + 1,2) = a;
        n = n + 1;
    end
    a = b;
    b = b - step_size;
    if abs(b) > abs((initial_guess - nmax)*step_size) % abs() since negative values
        break
    end
end
end

%% Bisection Method Function

% Bisection Pseudocode:
% Need to call the bracketing function to use the values that it obtains...
% since it's literally half of this method. Will need to define a lower and
% upper bound. They will be called a and b. The midpoint calculation will
% be called c. Will need an if statement to determine whether a or b should
% be redefined as c depending on the sign found. Then just have that repeat
% until the root is greater than the tolerance called for. 

function [roots_bisection, iteration_bisection, error_bisection] = bisection(root_guess, roots, f, tol)
roots_bisection = zeros(root_guess,1);

% a = x lower bound
% b = x upper bound
% c = x midpoint

for n = 1:root_guess
    a = roots(n,1);
    b = roots(n,2);
    iterations = zeros(root_guess,1);
    if sign(f(a))*sign(f(b)) >= 0
        error('f(a)f(b)<0 not satisfied!') % ceases execution
    end

    fa=f(a);
    while abs(b-a)/2>tol
        c=(a+b)/2; 
        fc=f(c);
        iterations(n,1) = iterations(n,1) + 1;
        q = iterations(n,1);
        iteration_bisection(n,q) = iterations(n,1);
        error_bisection(n,q) = c;
        if fc == 0 % c is a solution, done   
            break
        end
        if sign(fc)*sign(fa)<0 % a and c make the new interval
            b=c;
        else % c and b make the new interval
            a=c;
        end
    end
        roots_bisection(n)=(a+b)/2; % new midpoint is best estimate
end
end

%% Newton's Method Function

% Newton Pseudocode:
% Call the bracketing method. Need to include the first derivative of the
% f function as well. Need to define a lower and upper bound, say it is
% called a and b with a midpoint of xo. xO will simply be xo(2). While it 
% is not needed, starting from the midpoint makes the calculations fair 
% since the bisection, newton, and halley method will be compared. For 
% maximum efficiency, one would just use the lower bound. Just input the 
% formula for the newton method afterwards and have it in a while loop 
% until the tolerance is met.

function [roots_newton, iteration_newton, error_newton] = newton(root_guess, roots, f, g, tol) % need to call g since newton uses it
roots_newton = zeros(root_guess,1); % stores the roots found

% a = x lower bound
% b = x upper bound
% xo = x midpoint

xO = 0;
xoo = 0;

for n = 1:root_guess
    a = roots(n,1); % defines where to place the roots
    b = roots(n,2);
    iterations = zeros(root_guess,1);
    if sign(f(a))*sign(f(b)) >= 0
        error('f(a)f(b)<0 not satisfied!') % ceases execution
    end

    xo = (a+b)/2;
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

%% Halley's Method Function

% Halley Pseudocode:
% Literally the same as the newton method except now we are calling f, g,
% and h. xo is, again, being redefined and iterated until it reaches the
% desired tolerance. This redefining is handled by xO. Not sure of what 
% else to say other than we are using the bracking function's roots again.

function [roots_halley, iteration_halley, error_halley] = halley(root_guess, roots, f, g, h, tol) % need to call g & h since halley uses both
roots_halley = zeros(root_guess,1); % stores the roots found

% a = x lower bound
% b = x upper bound
% xo = x midpoint

xO = 0;
xoo = 0;

for n = 1:root_guess
    a = roots(n,1); % defines where to place the roots
    b = roots(n,2);
    iterations = zeros(root_guess,1);
    if sign(f(a))*sign(f(b)) >= 0
        error('f(a)f(b)<0 not satisfied!') % ceases execution 
    end

    xo = (a+b)/2;
    while abs(f(xoo))>tol % the halley formula
        xoo = xO; % whatever xO is, set equal to xoo for last iteration
        xO = xo - (2*f(xo)*g(xo))/(2*(g(xo))^2-f(xo)*h(xo));
        xo = xO; % redefine xo and xO to run the loop once more
        iterations(n,1) = iterations(n,1) + 1;
        e = iterations(n,1);
        iteration_halley(n,e) = iterations(n,1);
        error_halley(n,e) = xO;
    end

    roots_halley(n) = xO; % xO is redefined to be the roots 
        
end
end