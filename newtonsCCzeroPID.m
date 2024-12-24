%% newtonsCCzeroPID
% Created for ELEC 341, 2022 winter term
% Project Part 2, Q10
% Goal: PID controller complex zeros tuning to get the best phase margin possible
% Date: December 8, 2022
% Function Sytax
% INPUT:
%       Z_orig: A scalar complex number which will serve as the initial
%           search point (complex scalar)
%       GH: A open loop transfer function representing the dynamics except
%       for the zero block for a PID controlled feedback system (LTI)
% OUTPUT:
%       Zret: A scalar complex number which serves to create the best
%           available Phase Margin (complex scalar)
%       PMret: A scalar angle number which is the best phase margin this
%           iteration of newtonsCCzeroPID() was able to produce (scalar, deg)

function [Zret PMret] = newtonsCCzeroPID(Z_orig,GH)
% right off the bat, I want to supress the warning messages thrown by
% margin() failing to activate
warning('off','MATLAB:colon:nonIntegerIndex');
Z0 = struct('orig',[],'up',[],'down',[],'right',[],'left',[]);
PM = struct('orig',[],'up',[],'down',[],'right',[],'left',[]);

Z0.orig = Z_orig;
Dz = CCzero(Z0.orig);
if(~test_margin(GH*Dz))
    [Ggm Phm] = margin(GH);
    PM.orig = Phm;
else
    error('Origin Zero invalid and cannot activate margin()');
end
% now we must track in the cardinal directions. Moves in increments
% corresponding to 10^(-gRain)
gRain = 0; % defines search area adjustment quantity around target point.
% spiral is a variable to track how many loop iterations we've gone
% through. I'm not sure what to set the limit to, so you can definately
% change it to whatever you like. But according to my intuition, 1000 is a
% good benchmark.
spiral = 1;

while(gRain<15 && spiral<1000)
    % TESTING WITH DYNAMIC GRANULARITY (if the origin gets displaced, the
    % Grain will be decreased and a larger jump will be used. else if the
    % origin is static, the grain will be increased and smaller jumps used
    %disp(['Grain is at: ',num2str(gRain)]);
    [Z0.up PM.up] = searchTrace(Z0.orig,GH,PM.orig,gRain,'up');
    [Z0.down PM.down] = searchTrace(Z0.orig,GH,PM.orig,gRain,'down');
    [Z0.right PM.right] = searchTrace(Z0.orig,GH,PM.orig,gRain,'right');
    [Z0.left PM.left] = searchTrace(Z0.orig,GH,PM.orig,gRain,'left');

    % Now we move the origin with the stored results
    % From the creation, if
    % dir = 1, origin reached (centerpoint reached)
    % dir = 2, go up
    % dir = 3, go down
    % dir = 4, go right
    % dir = 5, go left
    PM_arr = [PM.orig PM.up PM.down PM.right PM.left];

    % Get the best phase margin which isn't infinity. That's probably from
    % an error or something in margin(). 
    [P_best dir] = max(PM_arr(~isinf(PM_arr)));
    if(dir~=1)
        %disp('Origin Reached, increasing granularity');
        gRain = gRain-1;
    else
        %disp('Origin Changed, decreasing granularity');
        gRain = gRain+1;
    end
    %P_best
    switch dir
        %case 1
            %disp('Origin Reached, increasing granularity');
            %gRain = gRain+1;
            % Case 1 is redundant when the above if() conditional will be
            % altering the gRain.
        case 2
            %disp(['Choose Zero: ',num2str(Z0.up)]);
            %disp('Go up');
            Z0.orig = Z0.up;
            PM.orig = PM.up;
        case 3
            %disp(['Choose Zero: ',num2str(Z0.down)]);
            %disp('Go down');
            Z0.orig = Z0.down;
            PM.orig = PM.down;            
        case 4
            %disp(['Choose Zero: ',num2str(Z0.right)]);
            %disp('Go right');
            Z0.orig = Z0.right;
            PM.orig = PM.right;            
        case 5
            %disp(['Choose Zero: ',num2str(Z0.left)]);
            %disp('Go left');
            Z0.orig = Z0.left;
            PM.orig = PM.left;            
    end
    Zret = Z0.orig;
    PMret = PM.orig;
    disp(['Best Zero: ',num2str(Zret)]);
    disp(['Best Phase Margin: ',num2str(PMret)]);
    % if the granularity gets too large, like 1e-15 being the prescision, you
    % can stop and just return whatever is in the Z0.orig
    spiral= spiral + 1;
end
%disp(['Spiral is at: ',num2str(spiral)]); % The quantity "spiral" tracks
%the number of loop iterations which have passed.
%disp("spiral limit reached");
warning('on','MATLAB:colon:nonIntegerIndex');
return;
end

%% Tracing Search for best PM improvement
% A function to search in any desired cardinal direction, signaled by a
% string. 'up','down','left','right'
% INPUTS:
%       Z0 to set as the origin (complex scalar)
%       GH open loop transfer function for PID without Dz block (LTI)
%       PM0 the origin's phase margin as a baseline to compare (real
%           scalar) in degrees
%       gRain is the degree of granularity. 1 means moving in steps of 1e-1
%       dir direction defining string (string)
% OUTPUTS:
%       bestZ0 greatest phase margin improvement zero in specified
%           direction within 2 of the original (complex scalar)
%       bestPM greatest phase margin (real scalar) in degrees

function [bestZ0 bestPM] = searchTrace(Z0,GH,PM0,gRain,dir)
bestZ0 = Z0; % Base zero
bestPM = PM0; % Base phase margin
stp = 0;
c=1;
while(stp<2)
    inc = c*10^(-gRain);
    switch dir
        case 'up'
            Z0_new = Z0+inc*1j;
        case 'down'
            Z0_new = Z0-inc*1j;
        case 'left'
            Z0_new = Z0-inc;
        case 'right'
            Z0_new = Z0+inc;
        otherwise
            error("invalid direction parameter");
    end
    Dz = CCzero(Z0_new);
    if(~test_margin(GH*Dz))
        stp = stp+1;
        [Ggm Phm] = margin(GH*Dz);
        if(Phm > bestPM)
            bestPM = Phm;
            bestZ0 = Z0_new;
            % a new zero was detected which yielded a better phase margin
        else
            return;
        end
    end
    if(c>50)
        error(['Error: margin has failed ', num2str(c),' times in a row.']);
    end
    c = c+1;
end
end
%% Test if margin() works
% internal flag. test_margin checks that the margin() function works. For
% some esoteric reason, margin() can fail randomly. It is essential to
% catch the error or else it will halt the program. smh.
% INPUT: sys (an open loop LTI system)
% OUTPUT: Boolean (0 or 1) if the margin() function can activate
%       1 means it can't
%       0 means it can
function C = test_margin(sys)
C = 0;
try
    bigboi = allmargin(sys);
catch E
    % Margin is unable to be initiated, this is a debugging piece of
    % code output
    %disp('Cannot initiate margin()');
    C = 1;
    % continue
end
end
%% Complex Conjugate Zero LTI block generator
% Take a complex zero and find the zeros block.
% INPUT: Z (complex zero)
% OUTPUT: Dzero (zeros block of a PID control system)
function Dzero = CCzero(Z)
s = tf('s');
Z_re = real(Z);
Z_im = imag(Z);
Z1 = Z_re + Z_im*1j;
Z2 = Z_re - Z_im*1j;
Dzero = 1/(Z1*Z2)*(s-Z1)*(s-Z2);
end
