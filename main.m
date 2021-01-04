%% Init
clear; clc;
fs = 100;
dt = 1/fs;
t = 0:dt:8-dt;
u =7*uGen(t, "step",1,100);
yraw = exciteSystem(5360188,u,fs);

% the laregest time constant was taken

%% Pre process

figure(1)
clf;
hold on; grid on;
plot(t,yraw);
plot(t,yrawl);
legend("ioanna","laurens");


%% Functions

function u = uGen(time,type, amp, periods)
	dt = time(2)-time(1);
	if type=="step"
		u = ones(length(time),1)*amp;
        u(1,1) = 0;
	elseif type == "pulse"
		u = [1 ; zeros(length(time)-1,1)];
	elseif type == "sine"
		u = amp*sin(periods*time*2*pi/(dt*length(time)));
	else 
		u = "Unknown input type";
	end
end