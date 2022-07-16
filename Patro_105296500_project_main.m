
% Author: Varun Patro
% Date: 08/08/20

% Clear Cache

clc; clear all; close all;

% defining variables
x_initial = 1;
v_initial = 0;
t_initial = 0;   
t_final = 10;
dt = 1/300;

% storing initial position and velocity
initial = [x_initial, v_initial];

% trial 1 conditions
m = 3;
k = 200;
c = 2;
omega = sqrt(k/m);
xi = (c)/(2*sqrt(m*k));
a = 2.5;

% setting up time scale for duration of timesteps
t = linspace(0,t_final,(ceil(t_final/dt))+1);

% defining variables for frequency response
xi_values = [sqrt(1/600), (45/2)*(sqrt(1/200)), 1, (35/2)*(sqrt(1/200)), 5*(sqrt(1/1000)),25*(sqrt(1/100)) ,30*(sqrt(1/900)), 2*(sqrt(1/900)) ];
n_exp = 8;
lambda = logspace(-1,1,500);
n_lambda = length(lambda);
Gain = zeros(n_exp, n_lambda);
PhiG = zeros(n_exp,n_lambda);
    
    % generating bode for all 8 trials
    figure(1)
    set(gcf,'Position',[15 50 1350 775])
    for run = 1:1:n_exp
        
        % implementing equations for gain and phase shift
        for j = 1:1:n_lambda
            Gain(run,j) = 1/sqrt((1-lambda(j)^2)^2 + (2*xi_values(run)*lambda(j))^2);
            PhiG(run,j) = -atan(2*xi_values(run)*lambda(j)/(1-lambda(j)^2));
        end
        % plotting with all information
        subplot(2,1,1)
            semilogx(lambda,Gain(run,:),'Linewidth',3)
            hold on
            grid on
            set(gca,'Linewidth',3,'Fontsize',20)
            ylabel('Amplification Factor')
            xlabel('Normalized Frequency')
            title('Function Gain')
          
        subplot(2,1,2) 
            semilogx(lambda,PhiG(run,:),'Linewidth',3)
            hold on
            grid on
            set(gca,'Linewidth',3,'FontSize',20)
            xlabel('Normalized Frequency')
            ylim([-2 2])
            ylabel('Radians')
            title('Function Phase Shift')
            
            legendInfo{run} = ['Trial ' num2str(run)];
    end
    legend(legendInfo);

n = input ('Enter 1 for Forward Euler, 2 for RK-2, or 4 for RK-4: ');

% each input corresponds to an approximation
if ((n == 1) || (n == 2) || (n == 4))
    switch n
        % Forward Euler Approximation
        case 1
            % setting up position arrays for inhomogenous and homogenous
            % response
            position = zeros(1,(ceil(t_final/dt))+1);
            v = VibrationPosition(initial,m,k,c,ForcingFunction(a,0),dt,1);
            % storing initial conditions inside
            position(1) = initial(1);
            position2 = zeros(1,(ceil(t_final/dt))+1);
            v2 = VibrationPosition(initial,m,k,c,0,dt,1);
            position2(1) = initial(1);
            % iterating through all time steps
            for b = 2:1:(ceil(t_final/dt))+1
                time = (b-1)*dt;
                f = ForcingFunction(a,time);
                f2 = 0;
                % calling function for each time step and setting equal to
                % itself
                v = VibrationPosition(v,m,k,c,ForcingFunction(a,time),dt,1);
                v2 = VibrationPosition(v2,m,k,c,0,dt,1);
                % storing values in arrays
                position(b) = v(1);
                position2(b) = v2(1);
                % plotting
                figure(2)
                plot(t,position,t,position2)
                xlabel('Time (s)')
                ylabel('Position (m)')
                title('Oscillatory Motion of Mass over Time (Forward Euler Model)')
            end
        case 2
            % contains same code above except for different approximation
            % is called
            position = zeros(1,(ceil(t_final/dt))+1);
            v = VibrationPosition(initial,m,k,c,ForcingFunction(a,0),dt,2);
            position(1) = initial(1);
            position2 = zeros(1,(ceil(t_final/dt))+1);
            v2 = VibrationPosition(initial,m,k,c,0,dt,2);
            position2(1) = initial(1);
            for b = 2:1:(ceil(t_final/dt))+1
                time = (b-1)*dt;
                f = ForcingFunction(a,time);
                f2 = 0;
                v = VibrationPosition(v,m,k,c,ForcingFunction(a,time),dt,2);
                v2 = VibrationPosition(v2,m,k,c,0,dt,2);
                position(b) = v(1);
                position2(b) = v2(1);
                figure(3)
                plot(t,position,t,position2)
                xlabel('Time (s)')
                ylabel('Position (m)')
                title('Oscillatory Motion of Mass over Time (RK 2 Model)')
            end
        case 4
            position = zeros(1,(ceil(t_final/dt))+1);
            v = VibrationPosition(initial,m,k,c,ForcingFunction(a,0),dt,4);
            position(1) = v(1);
            position2 = zeros(1,(ceil(t_final/dt))+1);
            v2 = VibrationPosition(initial,m,k,c,0,dt,4);
            position2(1) = v2(1);
            for b = 2:1:(ceil(t_final/dt))+1
                time = (b-1)*dt;
                f = ForcingFunction(a,time);
                f2 = 0;
                v = VibrationPosition(v,m,k,c,ForcingFunction(a,time),dt,4);
                v2 = VibrationPosition(v2,m,k,c,0,dt,4);
                position(b) = v(1);
                position2(b) = v2(1);
                figure(4)
                plot(t,position,t,position2)
                xlabel('Time (s)')
                ylabel('Position (m)')
                title('Oscillatory Motion of Mass over Time (RK 4 Model)')
            end 
    end
     
    % video file name
    video_filename = 'Underdamped';
    % setting up time scale
    timescale = t_initial:dt:t_final;
    
    % creating file
    vidfile = VideoWriter(video_filename,'MPEG-4');
    vidfile.FrameRate = 30;
    open(vidfile);
    % initializing frame
    nv = ceil(length(timescale)/10);
    t_out = zeros(nv,1);
    x_hom = zeros(nv,1);
    x_inh = zeros(nv,1);
    count = 0;
    sl = 0.27;
    
    figure(5)
    set(gcf,'Position',[15 50 1350 775])
    
    % iterating to make the frame update each timestep
    for i = 1:1:length(timescale)
        if mod(i-1,10) == 0
            count = count + 1;
            t_out(count) = t(i);
            x_hom(count) = position2(i);
            x_inh(count) = position(i);
            
            y1 = [x_hom(count)+sl x_hom(count)-sl x_hom(count)-sl x_hom(count)+sl];
            y2 = [x_inh(count)+sl x_inh(count)-sl x_inh(count)-sl x_inh(count)+sl];
            x = [-sl -sl sl sl];
            
            subplot(1,2,1)
                fill(x,y1,'b')
                xlim([-1-sl 1+sl])
                ylim([-1-sl 1+sl])
                xlabel('X location [m]')
                ylabel('Y location [m]')
                title('Homogenous Response')
                axis square
            subplot(1,2,2)
                fill(x,y2,'b')
                xlim([-1-sl 1+sl])
                ylim([-1-sl 1+sl])
                xlabel('X location [m]')
                title('Inhomogenous Response')
                axis square
                
                % saves video
                writeVideo(vidfile,getframe(gcf));
        end
    end
                
            
    
else 
    error('Invalid input. Enter 1,2, or 4');  
% End statement for very first if statement
end
        

