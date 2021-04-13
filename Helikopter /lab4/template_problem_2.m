% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Fl�ten

%% Initialization and model definition
init; % Change this to the init file corresponding to your helicopter

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A = [0 1 0 0 0 0;0 0 -K_2 0 0 0; 0 0 0 1 0 0;0 0 -K_1*K_pp -K_1*K_pd 0 0;0 0 0 0 0 1;0 0 0 0 -K_3*K_ep -K_3*K_ep];
I = eye(6);
A1 = I + delta_t*A;
%A1 = [1 delta_t 0 0 0 0;0 1 -K_2*delta_t 0 0 0; 0 0 1 1*delta_t 0 0;0 0 -K_1*K_pp*delta_t 1-K_1*K_pd*delta_t 0 0;0 0 0 0 1 delta_t;0 0 0 0 -delta_t*K_3*K_ep 1-delta_t*K_3*K_ep];
B1 = delta_t*[0 0;0 0;0 0;K_1*K_pp 0;0 0;0 K_3*K_4];

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               %e
x6_0 = 0;                               %e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';            % Initial values

% Time horizon and initialization
N  = 40;                               % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization
z0(1) = pi;

% Bounds
ul 	    = [-30*pi/180; -Inf];                            % Lower bound on control
uu 	    = [30*pi/180; Inf];                    % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                           % Lower bound on state x3
xu(3)   = uu(1);                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 2;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;
Q1(6,6) = 0;
q1 = 1;
q2 = 1;
R = [q1 0;0 q2];
P1 = 2*R;                                  % Weight on input
G = gen_q(Q1,P1,N,M);                    % Generate Q, hint: gen_q
%c = [];
% Generate c, this is the linear constant term in the QP

%% Generate system matrixes for linear model
Aeq = gen_aeq(A1,B1,N,mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(N*mx,1);                      % Generate b
beq(1:6)= A1*x0;
%% Solve QP problem with linear model
f = @(z) 1/2*z'*G*z;

tic
%[z,lambda] = quadprog(Q,c,[],[],Aeq,beq,vlb,vub,x0) ; % hint: quadprog. Type 'doc quadprog' for more info 
opt = optimoptions('fmincon', 'algorithm','sqp');
[z, fval] = fmincon(f,z0, [], [], Aeq, beq, vlb, vub, @constraints,opt);
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+G(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end


%% Lab 3
Q_lqr = diag([25 10 1 10 40 20]);
R_lqr = diag([1 1]);
K_lqr = dlqr(A1,B1,Q_lqr,R_lqr);


%% Extract control inputs and states
% u1  = [z(N*mx+1:2:N*mx+M*mu);z(N*mx+M*mu-1)]; % Control input from solution
% u2  = [z(N*mx+2:2:N*mx+M*mu);z(N*mx+M*mu)];
u1 = [z(N*mx+1:mu:N*mx+M*mu);z(N*mx+M*mu-1)];
u2 = [z(N*mx+2:mu:N*mx+M*mu);z(N*mx+M*mu)];

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution          

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];



%% Plotting
t = 0:delta_t:delta_t*(length(u1)-1);
% u1_opt = timeseries(u1,t); % Making time series for Simulink
% u2_opt = timeseries(u2,t); % Making time series for Simulink
% x1_opt = timeseries(x1,t); % Making time series for Simulink
% x2_opt = timeseries(x2,t); % Making time series for Simulink
% x3_opt = timeseries(x3,t); % Making time series for Simulink
% x4_opt = timeseries(x4,t); % Making time series for Simulink
% x5_opt = timeseries(x5,t); % Making time series for Simulink
% x6_opt = timeseries(x6,t); % Making time series for Simulink

x_opt = timeseries([x1,x2,x3,x4,x5,x6],t);
u_opt = timeseries([u1,u2],t);

figure(2)
sgtitle('Optimal trajectory, q = 1');


subplot(711)
stairs(t,u1),grid
ylabel('u')
subplot(712)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(713)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(714)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(715)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(716)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('elev')
subplot(717)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('elevdot')


function [c, ceq] = constraints(z)
   
    
    N = 40; mx = 6; alpha = 0.2; beta = 20; lambda_t = 2*pi/3;
    
    c = zeros(N,1);
    for k=1:N
        c(k) = alpha*exp(-beta*(z(1+(k-1)*mx)-lambda_t)^2)-z(5+(k-1)*mx);
    end
    ceq = [];
end