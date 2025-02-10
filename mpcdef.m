close all
clear
clc
%% Parametri del modello (dal paper)
m_s = 395.3;      % Massa sospesa (kg)
m_u = 48.3;       % Massa non sospesa (kg)
k_s = 30.01e+3;   % Rigidità sospensione (N/m)
c_s = 1450;     % Smorzamento (Ns/m)
k_t = 3.4e5;      % Rigidità pneumatico (N/m)

% Inputs: u = [F_act dot_zr]                   
% States: x = [zs-zu; dot_zs; zt; dot_zu]    
% Output: y= [a_s; susTravel; zt]
% y(1) = a_s = dotdot_zs = [-k_s/m_s -c_s/m_s k_s/m_s c_s/m_s]x + [1/m_s 0]u
% y(2) = susTravel = zs - zu = [1 0 -1 0]x +[0 0]u
% y(3) = zt = zu-zr = [0 0 1 0]x + [0 -1]

% Equation: dot_x = Ax + Bu
% Measurements: y = Cx + Du 

% A = [0, 1, 0, -1; -k_s/m_s, -c_s/m_s, 0, c_s/m_s; 
%      0, 0, 0, 1; k_s/m_u, c_s/m_u, -k_t/m_u, -c_s/m_u];
% B = [0, 0; 1/m_s, 0; 0, -1; -1/m_u, 0];
% C = [-k_s/m_s, -c_s/m_s, 0, c_s/m_s; 
%      1, 0, 0, 0; 0, 0, 1, 0];
% D = [1/m_s, 0; 0, 0; 0, 0];

A = [0 1 0 -1;-k_s/m_s -c_s/m_s 0 c_s/m_s;...
    0 0 0 1;k_s/m_u c_s/m_u -k_t/m_u -c_s/m_u];
B = [0 0;-1/m_s 0;0 -1;-1/m_u 0];    
C = [-k_s/m_s -c_s/m_s 0 c_s/m_s;...
    1 0 0 0;...
    0 0 1 0];                             
D = [1/m_s 0;0 0;0 0];
%% Discretizzazione del sistema
ts = 0.01; % Tempo di campionamento (s)
sys_d = c2d(ss(A, B, C, D), ts,'zoh');
[A_d, B_d,C_d,D_d] = ssdata(sys_d);

%% Parametri MPC
T_sim = 10;
T_c = 6;
T_p =12; % Orizzonte di previsione
% Q = diag([100, 2, 0.1, 2]);
Q = diag([0.1, 5, 0.1, 5]);
% Q = diag([0.1, 5, 0.1, 5]); % Peso sugli stati
% R = 1e-7; % Peso sugli ingressi
R = 0.0000001; % Peso sugli ingressi
[K, P] = dlqr(A_d, B_d,Q,R);
%calculates the optimal gain matrix K such that the state-feedback law
%minimizes the quadratic cost function; P is the infinite horizon solution 
%of the associated discrete-time Riccati equation
K = -K; % Sign convention.
Q_f= Q+K'*R*K;
% Q_f= Q;

x_lb = [-0.088;-0.163;-0.0128;-1.965]; %state lower bound
x_ub = [0.105;0.14;0.0128;2.78];       %state upper bound
u_lb = -2500; % Limiti inferiori sugli ingressi
u_ub = 2500; % Limiti superiori sugli ingressi

%% Profilo stradale (bump)
V = 30 * (1000 / 3600); % Velocità del veicolo (m/s)
time = 0:ts:T_sim; % Durata del profilo stradale (s)
bump_L = 5; % Lunghezza del bump (m)
bump_A = 0.1; % Altezza del bump (m)
t0 = 0.6; % Inizio del bump (s)

road_disturbance = zeros(1, length(time));
for t_idx = 1:length(time)
    t = time(t_idx);
    if t >= t0 && t <= t0 + bump_L / V
        road_disturbance(t_idx) = (bump_A / 2) * (1 - cos(2 * pi * V * (t - t0) / bump_L));
    end
end

%% Simulazione senza attuatore (passivo)
x_passive = zeros(4, length(time)); % Inizializzazione stati
y_passive = zeros(3, length(time)); % Inizializzazione uscite
u_control = zeros(1, length(time)); % Nessun controllo
u_disturbance = road_disturbance; % Disturbo stradale

for k = 1:length(time)-1
    % Dinamica del sistema discreto
    u_k = [u_control(k); u_disturbance(k)]; % Ingressi: controllo e disturbo
    x_passive(:, k+1) = A_d * x_passive(:, k) + B_d * u_k;
    y_passive(:,k) = C_d * x_passive(:,k) + D_d * u_k;
end
%% Simulazione MPC
T = length(time) - T_c;  % Orizzonte di simulazione
x = zeros(4, T+1); % Stati iniziali
x(:, 1) = [0; 0; 0; 0]; % Condizione iniziale
u_rec = zeros(T, 1); % Azioni di controllo (solo F_act)

for k = 1:T
    x_0 = x(:, k);  % Stato attuale
    z_road = road_disturbance(k:k+T_c-1)';  % Profilo stradale previsto (disturbo)

    %% Costruzione delle matrici di predizione T e S
    nx = 4;  % Numero di stati
    nu = 1;  % Ottimizzare solo F_act (un ingresso)
    T_mat = zeros(T_c*nx, nx);  % Matrice di evoluzione libera
    S_mat = zeros(T_c*nx, nu*T_c);  % Matrice di evoluzione forzata

    for i = 1:T_c
        T_mat((i-1)*nx+1:i*nx, :) = A_d^i;
        for j = 1:i
            S_mat((i-1)*nx+1:i*nx, j) = A_d^(i-j) * B_d(:, 1);  % Solo F_act
        end
    end

    %% Vincoli sugli stati in funzione di u_1

    % Vincoli: F * u <= e
    FS = [S_mat; -S_mat];  % Vincoli sugli ingressi
    e = [repmat(x_ub, T_c, 1) - T_mat * x_0;
         -repmat(x_lb, T_c, 1) + T_mat * x_0];

    %% Risoluzione del problema di ottimizzazione
    u0 = zeros(T_c, 1);  % Inizializzazione ingressi (solo F_act)
    lb = repmat(u_lb, T_c, 1);  % Limiti inferiori sugli ingressi
    ub = repmat(u_ub, T_c, 1);  % Limiti superiori sugli ingressi

    % Risoluzione con fmincon
    [u, fval] = fmincon(@(u_1) mpc_cost(u_1, A_d, B_d, Q, R, Q_f, x_0, z_road, T_c, T_p), ...
                            u0, FS, e, [], [], lb, ub, []);

    % Applica il primo ingresso ottimale
    u_rec(k) = u(1);  % Salva solo F_act (u_1)

    % Aggiorna stato
    x(:, k+1) = A_d * x_0 + B_d * [u_rec(k); z_road(1)];
    y(:,k) = C_d * x_0 + D_d * [u_rec(k); z_road(1)];
    
end
RMS_acc_pass = sqrt(1/100*(sum((y_passive(1,1:100)))^2))
RMS_acc_att = sqrt(1/100*(sum((y(1,1:100)))^2))
%% Risultati
time_mpc = time(1:T+1); % Tempo per simulazione con attuatore

figure;
for i = 1:nx
    subplot(nx, 1, i);
    hold on;
    plot(time_mpc, x(i, :), 'b', 'LineWidth', 1.5); % Stato attuale
    plot(time_mpc, repmat(x_lb(i), 1, length(time_mpc)), 'r--', 'LineWidth', 1); % Limite inferiore
    plot(time_mpc, repmat(x_ub(i), 1, length(time_mpc)), 'g--', 'LineWidth', 1); % Limite superiore
    hold off;

    ylabel(['x_' num2str(i)]);
    xlabel('Time (s)');
    title(['State x_' num2str(i) ' Evolution']);
    legend('State', 'Lower Bound', 'Upper Bound');
end


%% Risultati
% Creazione del vettore temporale corretto per la simulazione MPC
time_mpc = time(1:T+1); % Tempo per simulazione con attuatore
time_passive = time; % Tempo per simulazione passiva

% Grafici
figure;
subplot(2, 1, 1);
hold on;
plot(time_mpc, x(1, :), 'b', 'LineWidth', 1.5); % Sospensione con attuatore
plot(time_passive, x_passive(1, :), 'r--', 'LineWidth', 1.5); % Sospensione senza attuatore
% plot(time_passive,road_disturbance, 'k', 'LineWidth', 1.5);
hold off;
ylabel('Suspension Stroke (m)');
xlabel('Time (s)');
title('Suspension Stroke with and without Actuator');
legend('Active Suspension', 'Passive Suspension', 'Road Profile');

subplot(2, 1, 2);
plot(time(1:T), u_rec(1:T), 'r', 'LineWidth', 1.5); % Aggiunti marker
ylabel('Actuator Force (N)');
xlabel('Time (s)');
title('Actuator Force (Discrete)');
grid on;
%% Calcolo dell'accelerazione della massa sospesa (Sprung Mass)
acc_passive = y_passive(1,1:T);  % Accelerazione sospensione passiva
acc_active = y(1,1:T);  % Accelerazione sospensione attiva (uso i primi T valori di x)

%% Calcolo della deflessione del pneumatico (Tyre Deflection)
tyre_deflection_passive = x_passive(3, 1:T);  % Pneumatico passivo
tyre_deflection_active = x(3, 1:T);  % Pneumatico attivo (usando i primi T valori)

%% Plot dell'accelerazione della massa sospesa
figure;
hold on;
plot(time(1:T), acc_passive, 'r--', 'LineWidth', 1.5); % Accelerazione sospensione passiva
plot(time(1:T), acc_active, 'b', 'LineWidth', 1.5); % Accelerazione sospensione attiva
hold off;
ylabel('Acceleration of Sprung Mass (m/s^2)');
xlabel('Time (s)');
title('Acceleration of Sprung Mass: Active vs Passive Suspension');
legend('Passive Suspension', 'Active Suspension');
grid on;

%% Plot della deflessione del pneumatico
figure;
hold on;
plot(time(1:T), tyre_deflection_passive, 'r--', 'LineWidth', 1.5); % Deflessione pneumatico passiva
plot(time(1:T), tyre_deflection_active, 'b', 'LineWidth', 1.5); % Deflessione pneumatico attiva
hold off;
ylabel('Tyre Deflection (m)');
xlabel('Time (s)');
title('Tyre Deflection: Active vs Passive Suspension');
legend('Passive Suspension', 'Active Suspension');
grid on;

%% Funzione di costo
function f = mpc_cost(u,A_d, B_d, Q, R,Q_f,x_0,z_road,T_c,T_p)
    f = 0;
    x = x_0;
    % x_bar = [0.13; 0; 0.012; 0];
for ii = 1:T_c

    my_sum_i = [];

    for ll = 0:ii-1
        my_sum_i(:,ll+1) = A_d^(ll)*B_d*[u(ii-ll);z_road(ii-ll)];
    end

    x = A_d^(ii)*(x) + sum(my_sum_i,2);
    f(ii) = R*u(ii)^2 + (x)'*Q*(x);


end

x_Tc = x;

for ii = T_c+1:T_p-1

    x = A_d^(ii)*x_Tc;
    f(ii) = (x)'*Q*(x);

end


x_Tp = A_d*x;
f = sum(f)+(x_Tp)'*Q_f*(x_Tp);
end
