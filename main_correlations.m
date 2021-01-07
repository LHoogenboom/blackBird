%% Filtering & Identification
% PracticalAssignment 2
% Laurens Hoogenboom #4609638
% Ioanna Chanopoulou #5360188
clear; close all; clc;

%% Assignment 1

%%%Input Signal
% The considered system is a black box to us. We can input an arbitrarily
% chosen signal, but we don't know anyhting about sensor noise nor the
% dynamics of our system. First, we will use a simple step input to develop
% some preprocessing in assignment 1. In assignment 2, two completely new
% signals will be generated for identification and validation of the
% system. The signal used in Assignment 1 will not be used afterwards. 


%%% Sampling frequency
% The peak time is at .87 seconds and the signal starts rising after a
% delay of .45 seconds. This means that the rise time is about .42 seconds.
% An appropriate sample rate would be to have 8 or 9 samples in this time period. 
% So the sample time interval shoud be .42/9=.0467 seconds. i.e. a sampling 
% frequency of 22Hz when rounded up. In hindsight, an announcement was made
% that the signal generation process is correllated to the sampling
% frequency and that we were allowed to eyeball a good frequency. Double the
% found frequency (44Hz) gives a nice workable result. In addition, using a
% high frequency when we have to guess one, is a good choice since in thi
% way we can also prevent the aliasing phenomenon.

% the peak of the rise after a time delay was determined to be ca. .8
% seconds. This time was multiplied by 10 and taken as an appropriate
% duration for the simulation.
fs = 44; %Hz
dt = 1/fs;
t_gen = 0:dt:8-dt; % Need different t vector for u because U gives different length y.
u =7*uGen(t_gen, 'step',1,0); % unit step input for system analysis

%%% Persistency of excitation
% Persistently exciting input signal U (generated from unit step u)
U = genU(u); % works only for unit step input

%STUDENTID: the lowest-valued student-id of the students in the group.
%However, we are not using the lowest one because it doesn't give a nice
%response. So, we tried the other student number and it produced a
%managable result. 
STUDENTID =5360188; % 4609638 < 5360188
yraw = exciteSystem(STUDENTID, U, fs);
t = 0:dt/(countZeros(u)+1):(8)-dt/(countZeros(u)+1); % time vector for U

figure(1);
clf; hold on; grid on;
plot(t,yraw)
plot(t_gen,50*u)
legend('y_{raw}','Step*50')
xlabel('time [s]')
ylabel('signal magnitude')
title('Unprocessed (raw) measured data')

y = despike(yraw,10000,fs); % Desipike output

% The signal was despiked by means of, first, flatlining the spikes, and
% then intepolating between the beginning of the flatline and the first
% value after the flatline. The decision to first flatline the spike was
% made beucause this was just the first iteration of the despiking process. However,
% it is useful for determining the the width of the spike (it is not
% necessarily 1 extreme measurement) and thus it was kept. The
% interpolation performed was linear, using a function that we devveloped
% ourselves, called "interp".

figure(2)
clf; grid on;
plot(t,y)
xlabel('time [s]')
ylabel('signal magnitude')
title('Despiked Measured Data')

%%%Time delay - Time shifting
y = timeshift(y,500,fs);

cutoff = length(yraw)-length(y);
t_shifted = t(cutoff+1:end);
clear cutoff;

% The time shifting was performed looking at the slope of the output
% signal. This cuts off any part of the signal where the slope reaches 500,
% but only is the signal value at that point is equal or greater than 0.

% the right value for the slope criterion during despiking and timeshifting
% was found through trial and error. Nonetheless, it appears reliable for
% the input signal it was finetuned for.

figure(3)
clf; hold on;
plot(t_shifted,y)
xlabel('time [s]')
ylabel('signal magnitude')
title('Timeshifted Measurement Data')

% Note the initial respones time is not reset to 0 in the plot for visual effect, but it is
% when actual system identification is performed

y = DCoffset(y); % removes DC Offset

figure(4)
clf; hold on;
plot(t_shifted,y)
xlabel('time [s]')
ylabel('signal magnitude')
title('Pre-Processed Measured Data')


%%% Linearity Check

fprintf('%0s | %10s \n','Input gain','IO gain')
for i = 1:10 % determine IOgain for input gain 1 till 20
	ut =i*uGen(t_gen, 'step',1,9);

% Persistently exciting input signal
	Ut = genU(ut); % works only for unit step input
	yt = exciteSystem(STUDENTID, Ut, fs);
	iogain = IOgain(ut,yt,fs);
	fprintf('%-10.0i | %2.3f \n', ut(end), iogain)
end
clear ut Ut table;

% The table below shows IO gains for a variety of input gains. The IO gains
% seem constant by aproximation, which indicates linearity of the system.
% The system roughly behaves as: y*Igain=IOgain*u*Igain --> y=IOgain*u.E.I
% scaling the input gain scales the output gain with the same factor,
% resulting in a roughly constant IO gain of ca. 58.

%%% Generating training and validation data for Assignment 2

% Completely new data will be generated and pre-processed for the
% identification and validation of the system.

% Now, in order to do system identification we need training data and
% validation data. We will generate these using different types of signals.
% The training data will be generated using a Pseudo Random Binary Sequence (prbs). The
% validation data will be the data generated with a step function, like
% before. The DC offset was not removed for either output, because it
% appears this did not give a reliable model identification for us. But as
% can be seen before, removal itself does work.

fs = 100; %Hz
dur = 10;

t = 0:(1/fs):dur-(1/fs); % same time vector for generating both sets

u_v = uGen(t,'step',35,0); % validation data is step input mag=35
y_v = exciteSystem(STUDENTID,u_v,fs);

% Preprocess validation data
y_v = despike(y_v,80000,fs);
y_v = timeshift(y_v,4000,fs);
% y_v = DCoffset(y_v);

b = length(u_v)-length(y_v);
u_v = u_v(1:end-b); % shorten u_v for reasons explained in assignment 1


u_t = 35*prbs(length(t),0.5); % training data is PRBS mag=35
y_t = exciteSystem(5360188,u_t,fs);

y_t = despike(y_t,80000,fs);
y_t = timeshift(y_t,4000,fs);
% y_t = DCoffset(y_t);

a = length(u_t)-length(y_t);
u_t = u_t(1:end-a); % shorten input like before

% Shorten time to fit with either data set (for identification script to work)
t_t = t(1:length(u_t))';
t_v = t(1:length(u_v))';

clear a b;

% clear; % nicedata.mat was used to have cosistent data, for consistency and streamlining of the process
% load('nicedata.mat')

figure(12)
clf; hold on;
plot(t_t,y_t)
plot(t_v,y_v)
xlabel('time [s]')
ylabel('signal magnitude')
title('Training and Validation Output Data')
legend('training output','validation output')

%% Assignment 2
%Since we cannot use a white-noise sequence as an input signal, the
%identification methods that can be used are PI-MOESP and PO-MOESP. We will
%use PI-MOESP here, because we do not know if we are dealing with white
%measurement noise and because it gives a higher VAF value.
method = 'pi-moesp';

%To determine the order, we look for a significant gap between a set of dominant sin-
%gular values that correspond to the dynamics of the system and a set
%of small singular values due to the noise. This gap appears to vary every
%time, which is logical since the measurements have unpredictable spikes and
%are subject to random noise. However, a model order of n=5 seems to give the
%best overall result over multiple tries. It can sometimes be proven unreliable
%though.

n = 5;
s = 50;%We choose s to be n < s << N

% subspaceID runs into a problem, because the timeshifted output is no
% longer the same length as the input. Even though there is a delayed
% response, the output signal was originally the same length as the input
% signal. This means that the end of the input signal was partially not
% accounted for in the output. To solve this, the amount of measurements
% that were cut off at the beginning of the output will also be cut off at
% the end of the input, to make the prepocessed IO data more representative of
% eachother.

[A,B,C,D,x0,sv] = subspaceID(u_t,y_t,s,n,method); 

figure(6)
semilogy(sv,'x')
grid on; grid minor;
xlabel('time [s]')
ylabel('magnitude')
title('Singular Values of Identified SS system')


y_hat_t = simsystem(A,B,C,D,x0,u_t,fs,t_t);
%identification VAF
dy = y_t - y_hat_t; 
VAF_id = max(0, 1 - ((1/length(y_t))*sum(dy.^2))/((1/length(y_t))*sum(y_t.^2)))*100

% The VAF of the training output and simulated output with training input
% is high. It varies from time to time, but is consistently above 95%

% As can be seen below, the simulated model does a good job of tracking the
% measured output. However, this means that noise will be modeled too in
% the identification step. It is unclear to us at this point how much
% sensor noise there is and how much might be erratic behaviour of the
% system. So looking at the singular value plot we believe that, on
% average, a system order of 5 is good. What we notice is that during the
% rise of the input the system is tracked very well. It is on average best
% tracked with a system order of 5. We think that during this rise, the
% system response is most dominant, compared to noise. So we decided to
% tune identification for this.

figure(21)
clf; hold on;
plot(t_t,y_t)
plot(t_t,y_hat_t)
xlabel('time [s]')
ylabel('signal magnitude')
title('Training Data')
legend('Measured output','Simulated output')

%% Assignment 3

%%% Compare the VAFs

y_hat_val = simsystem(A,B,C,D,x0,u_v, fs,t_v);
%identification VAF
dy = y_v - y_hat_val; 
VAF_val = max(0, 1 - ((1/length(y_v))*sum(dy.^2))/((1/length(y_v))*sum(y_v.^2)))*100

% It appears as though the system still tracks well during the rise, but
% rapidly decays towards a steady state response. Because we use a
% completely different input this time the system might react differently.
% It is interesting that, though the simulated response in now very
% different, the VAF for the validation is still very hight, if not
% sometimes higher than the VAF for the training data. This might be
% coincedence, since a unit step and and prbs signal do resemble each other
% quite a lot.

figure(13)
clf; hold on;
plot(t_v,y_v)
plot(t_v,y_hat_val)
xlabel('time [s]')
ylabel('signal magnitude')
title('Validation Data')
legend('Measured output', 'Simulated output')

%%% Autocorrelation of residuals 
% We can evaluate model quality by checking whether the estimated auto-
% correlation function represents a unit pulse.
% epsilon: error-vector between model output and measured output for the 
% identification data set
epsilon1  =  y_t - y_hat_t;
epsilon2 = y_v - y_hat_val;
%lag interval d <= dur:experiment duration 
d  = 20;

%autocorrelation
ac1 = xcorr(epsilon1,epsilon1,d);
ac2 = xcorr(epsilon2,epsilon2,d);

%a plot of the autocorrelations follow
figure(15);plot(-d:d,ac1);hold on;plot(-d:d,ac2);hold off
xlabel('time lag')
legend('Autocorrelaion_{id}','Autocorrelaion_{val}')
title('Autocorrelation of the residuals')

%%%Cross-correlation between the residuals and the inputs 
% In case the estimated cross-correlation function equals zero, we have a
% good model estimation. As we can see from the following plot, the cross
% correlation we get is oscilating around the value of zero, so our
% model is sufficiently satisfying. 
 xc1 = xcorr(epsilon1,u_t,d);
 xc2 = xcorr(epsilon2,u_v,d);
%A plot can be made using the following commands
figure(11);plot(-d:d,xc1);
hold on;plot(-d:d,xc2);hold off
legend('Cross-correlation_{id}','Cross-correlation_{val}')
title('Cross correlation between the residuals and the inputs')


%% Functions used

function u = uGen(time,type, amp, periods) % generate inputs
	dt = time(2)-time(1);
	if type=="step" % gen step with periods as prior zeros
		u = [zeros(periods,1) ; ones(length(time)-periods,1)]*amp;
	elseif type == "pulse"
		u = [1 ; zeros(length(time)-1,1)];
    elseif type == "sine"
		u = amp*sin(periods*time*2*pi/(dt*length(time)));
	else 
		u = "Unknown input type";
	end
end

function u = prbs(N,rate)
    % wrote prbs because not allowed to use SI toolbox
    if rate>0 && rate<1
        du = floor(rand(N,1)+rate);
        u = rem(cumsum(du),2);
    else
        u = 'prbs(N,rate) | Error: Pick a rate between 0 and 1 \n';
        fprintf(u);
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
    rank = z+1; % for unit step input
    gain = max(u); % same input gain as 1D signal
    U = gain*ones(rank,length(u));
    U(1,:) = u'; % first row equals 1D signal
    for i=2:rank
		% Do bitshift operation (1 to the left) for all other rows with
		% respect to previous row.
       U(i,:) = [u(i:end)' gain*ones(1,i-1)]; 
    end
end

function y = interp(sig)
	% Interpolates flatlines signal until not flatlined anymore and appends
	% anything that comes after it. Start at flatline, will not influence
	% other flatlines
	done = false;
	i = 2;
	while ~done
		if ~(sig(i) == sig(1)) % if no longer flatlined
			dy = sig(i)-sig(1); % determine slope
			dx = i;
			for j=2:i-1 % for flatline
				sig(j) = sig(j-1)+dy/dx; % linear interpolation
			end
			done = true;
		end
		i = i+1; % look further if still at flatline
	end
	y = sig; % return
end

function y = despike(sig,slope,fs)
	% removes spikes by flatlining and then interpolating over flatline
	spikes = []; % to stare spike starts
    for i = 2:length(sig)
		
        if fs*(sig(i)-sig(i-1))>slope % check slope criterion
			
            sig(i) = sig(i-1); % flatline spike
			
			% only save start of identified spike
			if  ~(ismember(i-3,spikes)) && ~(ismember(i-2,spikes)) 
				spikes = [spikes i-1];
			end
        end
	end
	
	for i = 1:length(spikes) % interpolate over spikes and prepend prior data
		sig = [sig(1:spikes(i)-1) ; interp(sig(spikes(i):end))];
	end
    y =sig; % return
end

function y = timeshift(sig,slope,fs)
	% removes any data before certain slope is achieved at y>=0
    i = 2;
	done = false;
	while ~done % while no delayed respones identified
        if fs*(sig(i)-sig(i-1))>slope && sig(i-1)>=0 % if high and steep enough
            sig = sig(i-1:end); % remove prior data
			done = true; % don't run while loop again	
		end
		i = i+1;
	end
    y =sig; % return
end

function y = DCoffset(y)
	% removes constant offset to (ignores signal rise)
	done = false;
	i=2;
	while ~done
		% find first peak after rise
		if y(i) > y(i-1) && y(i) > y(i+1) && y(i) >400 
			sig = y(i:end); % look only after this peak
			done = true;
		end
		i = i+1;
	end
	DC = mean(sig); % remove mean after this peak from signal
	y = y-DC;
end

function g = IOgain(u,y,fs)
	% determines IO gain
	y = despike(y,10000,fs);
	y = timeshift(y,500,fs);
	g = mean(y)/mean(u);
end

function y_hat = simsystem(A,B,C,D,x0,u, fs,t)
%% Instructions:
% Implement a function that simulates the system here!
% Use the following function inputs and outputs.

% Function INPUT 
% A         System matrix A (matrix of size n x n)
% B         System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)
% u         system input (matrix of size N x m)

% Function OUTPUT
% yhat      predicted output (vector of size l x one)

%discrete time state space model, unknown sampling time (-1)
dt = 1/fs;
sys = ss(A,B,C,D,dt);

% %Estimation of the output
% Tfinal = (length(u)-1);
% t = 0:1/fs:Tfinal; %time

y_hat = lsim(sys,u,t,x0);

end

function [A,B,C,D,x0,sv] = subspaceID(u,y,s,n,method)
%% Instructions:
% Implement your subspace ID methods here.
% Avoid duplicate code! 
% Write the method specific code into the switch case!
% Use the following function inputs and outputs.

% Function INPUT 
% u         system input (matrix of size N x m)
% y         system output (matrix of size N x l)
% s         block size (scalar)
% n         model order (scalar)
% method    method (string e.g. 'moesp')
%
% Function OUTPUT
% A         System matrix A (matrix of size n x n)
% B         System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)
% sv        Singular values (vector of size n x one)
    m=1;%inputs
    l=1;%outputs
    switch method
        case 'moesp'
            %Step 1 - LQ factorization
            Han_u = hankel(u(1:l*s),u(s*l:end));%block Hankel matrix of the input 
            Y = hankel(y(1:l*s),y(s*l:end));%block Hankel matrix of the output sequence y
            L = (triu(qr([Han_u ; Y]')))';
            
            %Step 2 - SVD of L22-determinition of n,A,C
            L22 = L(s+1:2*s,s+1:2*s);
            [Un,Sn,Vn] = svd(L22);
            sv = diag(Sn);%Singular values
            A = Un(1:s*l-l,1:n)\Un(l+1:end,1:n);
            C = Un(1:l,1:n);
            
            %Step3 - Determine x0,B,D
            I = eye(l); 
            Phi = zeros(length(y),n+n*m +l*m);
            for k =0:(length(y)-1)
                SUM =zeros(1,n);
                for j =0:k-1
                    SUM = SUM + kron(u(j+1),C*A^(k-1-j));
                end
                Phi(k+1,:) = [C*A^k  SUM  kron(u(k+1,1),I)];
            end
            BD = Phi\y;
            x0 = BD(1:n);
            B = BD(n+1:2*n);
            D = BD(2*n+1);
        case 'pi-moesp'
            Up = hankel(u(1:m*s),u(s*m:length(y)-s));%past input 
            Uf = hankel(u(s+1:2*s), u(2*s*m:end)); %future input
            Yf = hankel(y(s+1:2*s),y(2*s:end)); %future output
            R =(triu(qr([Uf;Up;Yf]')))'; %LQ factorization
            R32 = R(2*s+1:2*s*m+s*l,s*m+1:s*m+s*l); %take R32 element
            
            [Un,Sn,Vn] = svd(R32);%SVD of R32
            sv = diag(Sn);%Singular values
            A = Un(1:s*l-l,1:n)\Un(l+1:end,1:n);
            C = Un(1:l,1:n);
            
            %Step3 - Determine x0,B,D
            I = eye(l); 
            Phi = zeros(length(y),n+n*m +l*m);
            for k =0:(length(y)-1)
                SUM =zeros(1,n);
                for j =0:k-1
                    SUM = SUM + kron(u(j+1),C*A^(k-1-j));
                end
                Phi(k+1,:) = [C*A^k  SUM  kron(u(k+1,1),I)];
            end
            BD = Phi\y;
            x0 = BD(1:n);
            B = BD(n+1:2*n);
            D = BD(2*n+1);
            
        case 'po-moesp'
            Up = hankel(u(1:m*s),u(s*m:length(y)-s));
            Uf = hankel(u(s+1:2*s), u(2*s*m:end));
            Yp = hankel(y(1:l*s),y(s*l:length(y)-s));
            Yf = hankel(y(s+1:2*s),y(2*s:end));  
            
            %-Step1 - LQ factorization
            Z =[Up;Yp];
            R =(triu(qr([Uf;Z;Yf]')))';
            R32 = R(3*s+1:3*s*m+s*l,s*m+1:s*m+s*l);% take R32
            
            %Step 2 - SVD
            [Un,Sn,Vn] = svd(R32);
            sv = diag(Sn);%Singular values
            A = Un(1:s*l-l,1:n)\Un(l+1:end,1:n);
            C = Un(1:l,1:n);
            
            %Step3 - Determine x0,B,D
            I = eye(l); 
            Phi = zeros(length(y),n+n*m +l*m);
            for k =0:(length(y)-1)
                SUM =zeros(1,n);
                for j =0:k-1
                    SUM = SUM + kron(u(j+1),C*A^(k-1-j));
                end
                Phi(k+1,:) = [C*A^k  SUM  kron(u(k+1,1),I)];
            end
            BD = Phi\y;
            x0 = BD(1:n);
            B = BD(n+1:2*n);
            D = BD(2*n+1);
    end
       
end