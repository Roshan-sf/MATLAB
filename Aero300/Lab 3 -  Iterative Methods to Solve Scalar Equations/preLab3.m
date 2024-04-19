%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 3 -  Iterative Methods to Solve Scalar Equations: 4/18/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Pseudo Code
% 
% rootMatrix = roots(initial_guess, step_size, expected_num_roots, function)
%     roots = [];  %create an empty matrix for later
%     current_point = initial_guess;
%     num_brackets = 0;
% 
%     while num_brackets < expected_num_roots
%         lower_bound = current_point;
%         upper_bound = current_point + step_size;
% 
%         %solve the function at lower and upper bounds
%         f_lower = function(lower_bound);
%         f_upper = function(upper_bound);
% 
%         % Check if the signs of bounds are different
%         if sign does not match
%             %found root interval [lower_bound, upper_bound]
%             roots = [roots; lower_bound, upper_bound];  %Append the bounds to the rootMatrix
%             num_brackets = num_brackets + 1;
% 
%             %Move the current point to the upper bound to find the next
%             bound
%             current_point = upper_bound;
%         else
%             %No bounds found, move to the next point
%             current_point = upper_bound;
%         end
%     end
% end
% 
% 
% 
%
