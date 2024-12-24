%% q21better

p=-CF/(N16);

Ts.old=Q20.Tr*.75;
OSu.old=Q20.OSu*.6;
inc=1;
c=0;

i=1*10^-inc;

tolerance=10^-inc;


GAIN.O.K = 1; 
GAIN.O.Kp = Q19.Kp; 
GAIN.O.Ki = Q19.Ki; 
GAIN.O.Kd = Q19.Kd; 

while abs(OSu.old-OSu.new)~=tolerance && abs(Ts.old-Ts.new)~=tolerance
    Q21.K=i;
    ip=1*10^-inc;
    while OSu.old-OSu.new~=0 && Ts.old-Ts.new~=0
        Q21.Kp=Q19.Kp*ip*Q21.K;
        ii=1*10^-inc;
        while OSu.old-OSu.new~=0 && Ts.old-Ts.new~=0
            Q21.Ki=Q19.Ki*ii*Q21.K;
            id=1*10^-inc;
            while OSu.old-OSu.new~=0 && Ts.old-Ts.new~=0
                Q21.Kd=Q19.Kd*id*Q21.K;
                GAIN.N.K = 1; % old 0.392
                GAIN.N.Kp = Q21.Kp; % old 0.1380
                GAIN.N.Ki = Q21.Ki; % old 0.3922
                GAIN.N.Kd = Q21.Kd; % old 0.0705

                [CLL OLL] = heurTune('PKID', 2, Q11.G, Q15.H, GAIN, p);

                info = stepinfo(CLL(5));
                peak = info.Peak;
                OSu.new = ((peak - 1) / 1) * 100;
                Ts.new = info.SettlingTime;
                if id > .1
                    break;
                end
                id=id+10^-inc;
            end
            if ii > .1
                break;
            end
            ii=ii+10^-inc;
        end
        if ip > .1
            break;
        end
        ip=ip+10^-inc;
    end
    if i>.1
        break;
        c=1;
    end
    i=i+10^-inc;
end

if c==1
    printf("no sol");
end
%disp(OSu.old-OSu.new);
%disp(Ts.old-Ts.new);

%Q21.Kp=GAIN.N.Kp;
%Q21.Ki=GAIN.N.Ki;
%Q21.Kd=GAIN.N.Kd;