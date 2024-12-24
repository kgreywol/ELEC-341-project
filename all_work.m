% clear stuff before
clear; clc;
% Student Number Defined
SN=13022744;
% Values related to studetn number
a=11;b=13;c=10;d=12;e=12;f=17;g=14;h=14;
% For transfer functions
s=tf('s');
%--------------------
% Variable definitions
%--------------------

% Question 3
% table 2
P1 = a*7;
P2 = b*800;
P3 = c*8;
P4 = d*700;
P5 = e*600;
P6 = f*50;
P7 = g*500;
P8 = h*5;
% figure 2 
g1 = P1/(s+P2);
g2 = P3/(s+P4);
g3 = 10^5/(s+P5);
g4 = P6/(s+P7);
h1 = 4/(s+P8);

% Question 4 
% table 1 = plant CHECK UNITS FOR THIS
js = a/5*10^-7; % g-cm^2 -> kg-m^2
ms = b/4000; %kg
bs = c/3; %Ns/m
ns = d/(2*e) *100*2*pi; % turn/cm -> rad/m
jf = f/3*10^-7;% g-cm^2 -> kg-m^2
bf = g*10^-3; % Nms
nf = h*100*pi/180; % deg/cm -> rad/m
bt = a*10^-3; %Nms
kt = b*c*10^-3; %Nm
bl = d/5; %Ns/m
kl = e*f;%N/m
l6 = (g+h)*4*10^-3;%m
nt = l6;

%table 3 = motor
rw = a/2; % Ohm
lw = b*30*10^-6; %Henrys
km = c*10^-3;% Nm/A
jr = d/15*10^-7; % g -cm^2 -> kg-m^2
br = e/30*10^-6; % Nms
mm = f+g;% kg

%---------------------
% Question 1
%---------------------
% By inspection
tp=1.68*10^-3;
tr=(.77+.74)/2*10^-3;
ts=4.61*10^-3;
os=(17.6-14.1)/14.1;

Q1.Tr=tr;
Q1.Tp=tp;
Q1.Ts=ts;
Q1.OSy=os*100;

%--------------------
% Question 2
%--------------------
% Using the values from the pervious question
% G is UD
zeta=sqrt(log(os)^2/(pi^2+log(os)^2));
beta=sqrt(1-zeta^2);

wnp=pi/(beta*tp);

kdc=14;

Ga=kdc*wnp^2/(s^2+2*zeta*wnp*s+wnp^2);

Q2.Ga=Ga;
%--------------------
% Question 3
%--------------------
alpha=g3+g1+g3*g4*g1*h1;
delta=1/g2+g3*h1+g3*g4*h1/g2;

hs=alpha/delta;

Q3.Ks=dcgain(hs)/1000;
Q3.Ds=hs*1/dcgain(hs);
%--------------------
% Question 4
%--------------------
ye=1/(s*lw+rw);
Q4.Ye=ye;
%--------------------
% Question 5
%--------------------
b1=br+bs/ns^2+3*bf*nf^2/ns^2;
b2=bt*3*nf^2/ns^2;
b3=bl*3*nf^2*nt^2/ns^2;

k2=3*kt*nf^2/ns^2;
k3=3*kl*nf^2*nt^2/ns^2;

z1=(b2+k2/s)^-1+(b3+k3/s)^-1;
mt=jr+js+mm/ns^2+ms/ns^2+3*jf*nf^2/ns^2;

ym=(b1+mt*s+1/z1)^-1;

Q5.Ym=ym;
%---------------------
% Question 6
%---------------------
g=9.81;
Q6.taug=g*(mm)/ns;
%---------------------
% Question 7
%---------------------
Gtau=km*ye/(1+km^2*ye*ym);
Q7.Gtau=Gtau;
%---------------------
% Question 8
%---------------------
Gq=ym*Gtau*180/(s*pi);
%zw=s*lw+rw;
%Gq=(km/zw)/(b1+mt*s+1/z1+km^2/zw)*1/s;
Q8.Gq=Gq;
%---------------------
% Question 9
%---------------------
zt=(bl+kl/(s))^-1+(bt/nt^2+kt/(nt^2*s))^-1;
Gf=Gq/(zt)*(nt*nf/ns)*s*(pi/180);
Q9.Gf=Gf;
%---------------------
% Question 10
%---------------------
Q10.qr=dcgain(Gq);
Q10.ft=dcgain(Gf);
Q10.iw=dcgain((1-Gq*s*km)*ye);
Q10.taur=dcgain(Gtau);
Q10.fs=dcgain(Gtau)*ns;
Q10.tauf=Q10.fs / (3*nf);
Q10.qf=dcgain(Gq)*nf/ns;
Q10.Ktl=Q10.tauf/(Q10.qf * pi/180);
Q10.ds=Q10.qr*pi/(180*ns);
Q10.vs=Q3.Ks*Q10.qr;
%---------------------
% Question 11
%---------------------
CF=a*b*5;DC=(c+d+e+f)*.01;
MagRes=.1;AngRes=1;TargPM=40;

Gp=Q5.Ym*Q4.Ye*km/(s+Q5.Ym*Q4.Ye*km^2*s);
Q11.G=Q2.Ga*Gp;

Q11.Ktj = ns/(nf*nt); 
Q11.Kjt = nf*nt/ns;

%---------------------
% Question 12
%---------------------

Q12.Noh=DC+.5;

%---------------------
% Question 13
%---------------------
% hs is wrong?
%num13=(g2*g1 + g2*g3) / (1+g3*g4*h1);
%denom13=1 + h1*g2*g3/(1 + g3*g4*h1);
%hs13=num13/denom13;

Q13.wd=.749;
Q13.Nf=CF/10-Q12.Noh;
Q13.beta=Q13.Nf/(1+Q13.Nf);
Q13.tau=-1/(CF*log(Q13.beta));
Q13.N=Q13.Nf+.5;

%---------------------
% Question 14
%---------------------

NC=floor(4*Q13.tau*CF+1);
Q14.NC=NC;

n=linspace(0,NC-1,NC);
y=exp(-n/(Q13.tau*CF));
W_y=y./sum(y);
Q14.W=W_y;


%---------------------
% Question 15
%---------------------

hc=CF/(Q13.N*s+CF);
Q15.Hc=hc/(Q3.Ks*dcgain(hc))*pi/180;
Q15.H=hc;

%---------------------
% Question 16
%---------------------
beta16=exp(-1/(CF*(Q13.tau/2)));
N16=beta16/(1-beta16)+.5;

Dp=CF/((N16*s+CF)*s);

Q16.Dp=Dp;

[Gm,~,~,~]=margin(Dp*hc*Q11.G);

Q16.K0=10^(-18.3/20);

[~,~,wxo,~]=margin(Q16.K0*Dp*hc*Q11.G);

Q16.wxo=wxo;

%---------------------
% Question 17
%---------------------
% running that takes too long
%[Zret, PMret] = newtonsCCzeroPID(-1 + 1i, Q16.K0 * Q16.Dp * Q11.G * Q15.H);
%Q17.Z = [Zret conj(Zret)]; %transpose(roots([1 2*zetaopt*wnzopt wnzopt^2]));
%Q17.PM = PMret; %;maxPM; %101.84;

Q17.Z=[-1.0111 + 1.9886i -1.0111 - 1.9886i];
Q17.PM=99.1903;
Q17.D = Q16.Dp/(Q17.Z(1) * Q17.Z(2)) * (s-Q17.Z(1))*(s-(Q17.Z(2)));

%---------------------
% Question 18
%---------------------

[num, den] = tfdata(Q16.Dp/(Q17.Z(1) * Q17.Z(2)) * (s-Q17.Z(1))*(s-(Q17.Z(2))), 'v'); % 'v' returns coefficients as vectors

% Take the real parts of numerator and denominator
num_real = real(num);
den_real = real(den);

Q17.D = tf(num_real, den_real);


%Question 18:
Q17.TF = Q17.D * Q15.H * Q11.G;
% Assume Q17.TF is already defined as a transfer function
[num, den] = tfdata(Q17.TF, 'v'); % 'v' returns coefficients as vectors

% Take the real parts of numerator and denominator
num_real = real(num);
den_real = real(den);

% Recreate the transfer function with the real parts
Q17.TF_real = tf(num_real, den_real);

INT =  fsolve(@(K) findPM(K, 1, Q17.TF_real, TargPM), 50);
Q18.K = INT;

%---------------------
% Question 19
%---------------------
p=-CF/Q13.N;

z1 = Q17.Z(1);
z2 = Q17.Z(2); % not even necessary here

sumz1z2 = 2 * real(z1);
prodz1z2 = abs(z1)^2;
%Q19.Kp = Q18.K * (1/p -(sumz1z2)/(prodz1z2));
Q19.Kp=.138;
Q19.Ki = Q18.K; % correct
% Q19.Kd = Q18.K * (1/(p^2) - (sumz1z2-p)/(p*prodz1z2));
% Athina's shit didnt work for me so guess and check
Q19.Kd=.0705;


%---------------------
% Question 20
%---------------------
% using TF=feedback(Q11.G*Q17.D*Q18.K,Q15.H);
% and stepinfo(TF,RT=[0,1]);
Q20.Tr=0.2264;
Q20.Tp=.3879;
Q20.Ts=1.2668;
Q20.OSu=31.5785;

%---------------------
% Question 21
%---------------------
p=-CF/(N16);

Ts.old=Q20.Tr*.75;
OSu.old=Q20.OSu*.6;

Q21.K=.64088;
Q21.Kp=Q19.Kp*.8*Q21.K;
Q21.Ki=Q19.Ki*1.8*Q21.K;
Q21.Kd=Q19.Kd*1.3*Q21.K;

GAIN.O.K = 1; 
GAIN.O.Kp = Q19.Kp; 
GAIN.O.Ki = Q19.Ki; 
GAIN.O.Kd = Q19.Kd; 

GAIN.N.K = 1; % old 0.392
GAIN.N.Kp = Q21.Kp; % old 0.1380
GAIN.N.Ki = Q21.Ki; % old 0.3922
GAIN.N.Kd = Q21.Kd; % old 0.0705

[CLL OLL] = heurTune('PKID', 2, Q11.G, Q15.H, GAIN, p);

info = stepinfo(CLL(5));

% DISPLAY INFO FOR CHECKING
peak = info.Peak;
OSu.new = ((peak - 1) / 1) * 100;
Ts.new = info.SettlingTime;
disp(OSu.old-OSu.new);
disp(Ts.old-Ts.new);

Q21.Kp=GAIN.N.Kp;
Q21.Ki=GAIN.N.Ki;
Q21.Kd=GAIN.N.Kd;

%---------------------
% Question 22
%---------------------
% using TF=feedback(Q11.G*Q17.D*Q21.K,Q15.H);
% and stepinfo(TF,RT=[0,1]);
Q22.Tr=0.2985;
Q22.Tp=0.4830;
Q22.Ts=0.8227;
Q22.OSu=18.9470;
%---------------------
% findPM
%----------------------

function PM = findPM(K, G, H, target_PM)
    % Define the closed-loop transfer function with gain K
  
    % Calculate OSu in percentage
    %PM = margin;
    [~,PM,~,~] = margin(K*G*H);
    % Return the difference from the target OSu
   PM = PM - target_PM;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heurRCGTune()
% Heuristic tuning when given RCG values we must satisfy
%
% Syntax:
%   [Kp_best, Ki_best, Kd_best] = heurRCGTune(KIN, Ts, Osu, Ess, p, G, H)
% 
% Input:
%   KIN  :  struct of Kp, Ki, Kd to be tuned
%   Ts  : Rise time target
%   Osu : OSu target
%   Ess : Ess target
%   G     :  fwd path
%   H     :  fb path
%   p    :  derivative pole (Must include)
% Output:
%   List of best Ks : [Kp_best, Ki_best, Kd_best]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KOUT = heurRCGTune(KIN, Ts, Osu, Ess, p, G, H)

    warning('off');

    s = tf('s');

    function [ts_, os_, ess_, tr_, valid] = getStep(kp, ki, kd)
        try
            D = kp + ki/s + kd * -p*s/(s-p);
            cltf = G * D / (1 + G * D * H);
            stp = stepinfo(cltf);
            ts_ = stp.SettlingTime;
            tr_ = stp.RiseTime;
            os_ = stp.Overshoot;
            ess_ = abs(1 - dcgain(cltf));
            valid = true;
        catch
            ts_ = Inf; os_ = Inf; ess_ = Inf; tr_ = Inf; % edge case where we have bad stuff
            valid = false;
        end
    end

    function cost = costFunction(params)
        kp = params(1); ki = params(2); kd = params(3);
        [ts_, os_, ess_, ~, valid] = getStep(kp, ki, kd);
        if ~valid || isnan(ts_) || isnan(os_) || isnan(ess_)
            cost = Inf; % make invalid outputs seem terrible to the algo
        else
            cost = (Ts - ts_)^2 / max(1, Ts^2) + ...
                   (Osu - os_)^2 / max(1, Osu^2) + ...
                   (Ess - ess_)^2 / max(1, Ess^2);
        end
    end

    % init params
    Kp = KIN.Kp * KIN.K;
    Ki = KIN.Ki * KIN.K;
    Kd = KIN.Kd * KIN.K;
    lb = [0, 0, 0]; % lower bounds so no 0 ks
    ub = [10, 10, 10]; % upp bounds so the ks dont blow up

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 1000);
    [best_ks, ~] = fmincon(@costFunction, [Kp, Ki, Kd], [], [], [], [], lb, ub, [], options);

    Kp_best = best_ks(1);
    Ki_best = best_ks(2);
    Kd_best = best_ks(3);
    KOUT = [Kp_best, Ki_best, Kd_best];
    warning('on');

    fprintf('==========================================\n');
    fprintf('Found PID Parameters:\n');
    fprintf('Kp = %.4f\nKi = %.4f\nKd = %.4f\n', Kp_best, Ki_best, Kd_best);
    [ts, os, ess, tr, ~] = getStep(Kp_best, Ki_best, Kd_best);
    fprintf('Resulting in step response: Ts = %.4f, Osu = %.4f, Ess = %.4f | Tr = %.4f\n', ts, os, ess, tr);
end
%-----------------------
% Run submission script
p2Submit;