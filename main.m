%% Init
clear; clc;
t = 0:.1:99;
u = 50*uGen(t, "step",1,100);
fs = 100;
yraw = exciteSystem(4609638,u,fs);

%% Pre process

figure(1)
clf;

y = lowpass(yraw,.001,fs);

figure(2)
clf;
hold on; grid on;
plot(t,u);
plot(t,y);
legend("input","output");


%% Functions

function u = uGen(time,type, amp, periods)
	dt = time(2)-time(1);
	if type=="step"
		u = ones(length(time),1)*amp;
	elseif type == "pulse"
		u = [1 ; zeros(length(time)-1,1)];
	elseif type == "sine"
		u = amp*sin(periods*time*2*pi/(dt*length(time)));
	else 
		u = "Unknown input type";
	end
end