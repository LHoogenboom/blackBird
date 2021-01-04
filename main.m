%% Init
clear; clc;
fs = 100; % Hz
dt = 1/fs;
t = 0:dt:8-dt;
u =7*uGen(t, "step",1,0);
yraw = exciteSystem(5360188,u,fs);

% The peak time is at .87 seconds and the signal starts rising after a
% delay of .45 seconds. This means that the rise time is about .42 seconds.
% An appropriate sample rate would be to have 8 or 9 samples in this time period. 
% So the sample time interval shoud be .42/9=.0467 seconds. i.e. a sampling 
% frequency of 22Hz when rounded up

% the peak of the rise after a time delay was determined to be ca. .8
% seconds. This time was multiplied by 10 and taken as an appropriate
% duation for the simulation.

figure(1);
clf; hold on; grid on;
plot(t,u*50)
plot(t,yraw)
legend("50 * Input_{step}","y_{raw}")

%% Peristence of excitation

m = 1; %number of inputs
l = 1; %number of outputs
sysOrder = 10;
U = genU(u); % works only for unit step input

%% spikes
y = despike(yraw,100);


figure(2)
clf; grid on;
plot(t,y)


%% Timeshift


%% Functions

function u = uGen(time,type, amp, periods)
	dt = time(2)-time(1);
	if type=="step"
		u = [zeros(periods,1) ; ones(length(time)-periods,1)]*amp;
	elseif type == "pulse"
		u = [1 ; zeros(length(time)-1,1)];
	elseif type == "sine"
		u = amp*sin(periods*time*2*pi/(dt*length(time)));
	else 
		u = "Unknown input type";
	end
end

function z = countZeros(u) % counts zeros before unit step
    z = 0;
    while u(z+1)==0
        z=z+1;
    end
end

function U = genU(u)
    % Generates Persistently exciting unit step 
    z = countZeros(u);
    sysOrder = z+1;
    gain = max(u);
    U = gain*ones(sysOrder,length(u));
    U(1,:) = u';
    for i=2:sysOrder
       U(i,:) = [u(i:end)' gain*ones(1,i-1)]; 
    end
end

function y = despike(sig,steep)
    for i = 2:length(sig)
        if sig(i)-sig(i-1)>steep
            sig(i) = sig(i-1);
        end
    end
    y =sig;
end

function y = timeshift(u,y)
    delay = .45;
    i = 1;
    while u(i)==0
        i=i+i;
    end
    y = y(i:end)
end
