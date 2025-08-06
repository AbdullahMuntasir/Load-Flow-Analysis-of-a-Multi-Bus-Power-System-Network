clc; clear; close all;

% Base values
S_base = 100;   % MVA
V_base = 138;   % kV

%% --- BUS DATA: [BusNo Type Pg(MW) Qg(MVAr) Pl(MW) Ql(MVAr) Vset(pu) Angle(deg)] ---
% Type: 1=Slack, 2=PV, 3=PQ

bus_data = [
    1 1   0    0      0    0    1.00  0;   % Slack (BUS 1)
    2 2  330  165     0    0    1.00  0;   % PV/Gen (BUS 2)
    3 2  330  165     0    0    1.00  0;   % PV/Gen (BUS 3)
    4 3   0    0    115   10    1.00  0;   % PQ/Load (BUS 4)
    5 3   0    0    115   10    1.00  0;   % PQ/Load (BUS 5)
    6 3   0    0    114   10    1.00  0;   % PQ/Load (BUS 6)
    7 3   0    0    115   10    1.00  0;   % PQ/Load (BUS 7)
    8 3   0    0    115   10    1.00  0;   % PQ/Load (BUS 8)
    9 3   0    0    115   10    1.00  0;   % PQ/Load (BUS 9)
];

n_bus = size(bus_data,1);

%% Convert MW/MVAr to per unit
bus_data(:,3:6) = bus_data(:,3:6) / S_base;

type = bus_data(:,2);
P = bus_data(:,3) - bus_data(:,5);
Q = bus_data(:,4) - bus_data(:,6);

%% --- LINE DATA: [From To R X] ---
line_data = [
    1 2 0.02 0.1;
    1 4 0.02 0.1;
    2 3 0.02 0.1;
    3 9 0.02 0.1;
    4 5 0.02 0.1;
    5 6 0.02 0.1;
    6 7 0.02 0.1;
    6 8 0.02 0.1;
    7 2 0.02 0.1;
    9 8 0.02 0.1;
];

%% --- Build Ybus ---
Ybus = zeros(n_bus, n_bus);
for k = 1:size(line_data,1)
    i = line_data(k,1);
    j = line_data(k,2);
    z = line_data(k,3) + 1j*line_data(k,4);
    y = 1/z;
    Ybus(i,i) = Ybus(i,i) + y;
    Ybus(j,j) = Ybus(j,j) + y;
    Ybus(i,j) = Ybus(i,j) - y;
    Ybus(j,i) = Ybus(j,i) - y;
end

%% --- Initial Values ---
V = bus_data(:,7) .* exp(1j * deg2rad(bus_data(:,8)));
max_iter = 1000;
tol = 1e-6;
iter = 0;

slack = find(type==1);

% PV bus Q-limits in p.u.
Qmax = 165 / S_base; % 0.2 p.u. (for this ID)
Qmin = -165 / S_base;

%% --- Gauss-Seidel Iteration with Q-limits ---
while iter < max_iter
    V_prev = V;
    iter = iter + 1;
    for i = 1:n_bus
        if type(i) == 1  % Slack
            continue;
        end
        sumYV = Ybus(i,:) * V - Ybus(i,i) * V(i);
        if type(i) == 3  % PQ
            S = P(i) + 1j*Q(i);
            V(i) = (1/Ybus(i,i)) * ((conj(S)/conj(V(i))) - sumYV);
        elseif type(i) == 2  % PV with Q limit
            Q_temp = -imag(V(i) * conj(Ybus(i,:) * V));
            if Q_temp > Qmax
                Q(i) = Qmax;
                type(i) = 3; % Convert to PQ
            elseif Q_temp < Qmin
                Q(i) = Qmin;
                type(i) = 3;
            else
                Q(i) = Q_temp;
                S = P(i) + 1j*Q(i);
                V(i) = (1/Ybus(i,i)) * ((conj(S)/conj(V(i))) - sumYV);
                V(i) = abs(bus_data(i,7)) * exp(1j * angle(V(i))); % fix |V|
            end
        end
    end
    if max(abs(V - V_prev)) < tol
        break;
    end
end

%% --- Output: Bus Voltages and Angles ---
fprintf('\nBUS VOLTAGE RESULTS AFTER %d ITERATIONS:\n', iter);
fprintf('Bus\t|V| (p.u.)\tAngle (deg)\t|V| (kV)\n');
for i = 1:n_bus
    Vmag_pu = abs(V(i));
    Vang_deg = rad2deg(angle(V(i)));
    Vmag_kV = Vmag_pu * V_base;
    fprintf('%d\t%.4f\t\t%7.2f\t\t%.2f\n', i, Vmag_pu, Vang_deg, Vmag_kV);
end
%% === Line Flow Calculation ===
branch = line_data(:,1:2);
num_lines = size(branch,1);
line_flows = zeros(num_lines, 1);

for k = 1:num_lines
    i = branch(k,1);
    j = branch(k,2);
    z = line_data(k,3) + 1j*line_data(k,4);
    y = 1/z;
    Iij = (V(i) - V(j)) * y;
    Sij = V(i) * conj(Iij) * S_base;
    line_flows(k) = real(Sij);  % Active power flow in MW
end

%% === Graphical Plots ===
figure('Name','Load Flow Results - ID 231033');

% --- Bus Voltage Magnitudes ---
subplot(3,1,1);
bar(abs(V), 'FaceColor', [0.1 0.7 0.4]);
title('Bus Voltage Magnitudes');
xlabel('Bus Number');
ylabel('|V| (p.u.)');
ylim([0 1.1]);
grid on;

% --- Voltage Angles ---
subplot(3,1,2);
stem(rad2deg(angle(V)), 'filled', 'Color', [0.8 0.2 0.2]);
title('Voltage Angles');
xlabel('Bus Number');
ylabel('Angle (°)');
grid on;

% --- Line Active Power Flows ---
subplot(3,1,3);
plot(line_flows, '-s', 'LineWidth', 2, 'Color', [0.3 0.3 1]);
title('Line Active Power Flows');
xlabel('Line Index');
ylabel('P (MW)');
grid on;
%% --- CALCULATE LINE FLOWS, LOSSES, AND SYSTEM LOSS ---
fprintf('\nLINE FLOWS (Sending-End):\n');
fprintf('From-To\tP_ij (MW)\tQ_ij (MVAr)\n');
total_loss = 0;

for k = 1:size(line_data,1)
    i = line_data(k,1);
    j = line_data(k,2);
    z = line_data(k,3) + 1j*line_data(k,4);
    y = 1/z;
    I_ij = (V(i) - V(j)) * y;
    S_ij = V(i) * conj(I_ij);  % Sending-end power (from i to j)
    I_ji = (V(j) - V(i)) * y;
    S_ji = V(j) * conj(I_ji);  % Sending-end power (from j to i)
    loss = S_ij + S_ji;

    fprintf('%d-%d\t%8.2f\t%8.2f\n', i, j, real(S_ij)*S_base, imag(S_ij)*S_base);

    total_loss = total_loss + loss;
end

%% --- PRINT LINE LOSSES ---
fprintf('\nLINE LOSSES:\n');
fprintf('From-To\tLoss (MW)\tLoss (MVAr)\n');
for k = 1:size(line_data,1)
    i = line_data(k,1);
    j = line_data(k,2);
    z = line_data(k,3) + 1j*line_data(k,4);
    y = 1/z;
    I_ij = (V(i) - V(j)) * y;
    S_ij = V(i) * conj(I_ij);
    I_ji = (V(j) - V(i)) * y;
    S_ji = V(j) * conj(I_ji);
    line_loss = S_ij + S_ji;
    fprintf('%d-%d\t%8.2f\t%8.2f\n', i, j, real(line_loss)*S_base, imag(line_loss)*S_base);
end

%% --- PRINT TOTAL SYSTEM LOSS ---
fprintf('\nTOTAL SYSTEM LOSS:\n');
fprintf('Total Loss (MW):   %.2f\n', real(total_loss)*S_base);
fprintf('Total Loss (MVAr): %.2f\n', imag(total_loss)*S_base);