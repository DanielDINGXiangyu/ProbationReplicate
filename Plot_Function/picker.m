function [indx] = picker(vect)
% This function pick up the state index of
%    Input: vect = [vect1 vect2 ...]
%       vect1 is a intiger, first digit is zlevel (1 lowest, 5 highest)
%       second digit is zvol, third digit is glevel, forth digit is gvol
%    Output: indx is the index of column of the state in statespace

    NZlevel=5;
    NZsigma=3;
    NGlevel=5;
    NGsigma=3;
    
    % [zlevel zvol glevel gvol]: output is the index of the state
    picker_scalar = @(x)...
        +NGsigma*NGlevel*NZsigma*(fix(mod(x,10000)/1000)-1)... % tlevel
        +NGsigma*NGlevel*(fix(mod(x,1000)/100)-1)...           % tvol
        +NGsigma*(fix(mod(x,100)/10)-1)...                     % glevel
        +(mod(x,10));
    
    % vectorized the function
    picker_vector =  @(vect) arrayfun(picker_scalar, vect);
    
    indx =  picker_vector(vect);
    
end

