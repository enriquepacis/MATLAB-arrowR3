close all;
clear all;

%% Constants

% Pauli matrices
sx = [0 1; 1 0];
sy = 1i*[0 -1; 1 0];
sz = [1 0; 0 -1];

hbar = 6.582E-16; % [eV*s, energy * time] the reduced Planck constant

ket0 = [1; 0];
ket1 = [0; 1];

%% Parameters
gamma = 0.05; % [eV]
Delta = 0;    % [eV] detuning

nt = 60; % [] number of time steps
T0 = pi*hbar/gamma; % [s] a convenient time scale
t_max = 2*T0; % [s] maximum calculation time

% Initial state
psi_0 = ket0;

arrowstyle = 'arrowR3()'; % here, use either 'arrowR3()' or 'arrow3()'


%% Calculations
t = linspace(0, t_max, nt); % [s] time vector

H = -gamma * sx + Delta/2 * sz;

Gx = trace(H * sx) / hbar; % [Hz]
Gy = trace(H * sy) / hbar; % [Hz]
Gz = trace(H * sz) / hbar; % [Hz]
Gamma = [Gx; Gy; Gz];

% Make a unit vector in the direction of Gamma

GammaDir = Gamma * (1/sqrt(Gamma' * Gamma));

psi_t = zeros(2, nt);
lx = zeros(1, nt);
ly = zeros(1, nt);
lz = zeros(1, nt);

for t_idx = 1:nt
    Ut = expm( -1i * H * t(t_idx)/hbar );
    psi_t(:, t_idx) = Ut * psi_0;
    
    rho_t =  psi_t(:, t_idx) * (psi_t(:, t_idx))';
    lx(t_idx) = real(trace(rho_t * sx));
    ly(t_idx) = real(trace(rho_t * sy));
    lz(t_idx) = real(trace(rho_t * sz));
end


%% Visualization

plot( t/T0, lz, 'LineWidth', 2)
grid on;
set(gca, 'FontSize', 18, 'FontName', 'Times')
xlabel('$t/T_0$', 'Interpreter', 'latex')
ylabel('$\lambda_z$', 'Interpreter', 'latex')

% Movie
figure;

M(nt) = struct('cdata', [], 'colormap', []);

for t_idx = 1:nt
    % Plot the Bloch sphere
    [Xs, Yx, Zx] = sphere(25);
    mySphere = surf( Xs, Yx, Zx );
    axis equal
    shading interp
    mySphere.FaceAlpha = 0.25
    set(gca, 'FontSize', 18, 'FontName', 'Times')
    
    line([-1 1], [0 0], [0 0], 'LineWidth', 1, 'Color', [0 0 0]); % x-axis
    line([0 0], [-1 1], [0 0], 'LineWidth', 1, 'Color', [0 0 0]); % y-axis
    line([0 0], [0 0], [-1 1], 'LineWidth', 1, 'Color', [0 0 0]); % z-axis
    
    O = [0 0 0]; % origin
    L = [lx(t_idx), ly(t_idx), lz(t_idx)]; % tip of coherence vector
    
    hold on
    plot3(lx, ly, lz, 'LineWidth', 2); % plot the vector's path
    switch arrowstyle
        case 'arrow3()'
            % arrow3(P1,P2,S,W,H,IP,ALPHA,BETA)
            arrow3(O, L, 'b', 3, 5, 4);
            arrow3(O, GammaDir', 'r', 3, 5, 4);
        case 'arrowR3()'
            arrowR3(O, L)
            arrowR3(O, GammaDir, 'ArrowHeadColor', [1, 0.6, 0.6], ...
                'ArrowHeadColor', 0.6*[0.8, 0.2, 0.2])
    end

    hold off
    view([-65, 10])
    title(arrowstyle)
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    zlabel('$z$', 'Interpreter', 'latex')
    drawnow;
    
    M(t_idx) = getframe(gcf);
    
    close(gcf); 
end

writerObj = VideoWriter([arrowstyle(1:end-2), '.mp4'], 'MPEG-4');
writerObj.FrameRate = 12;

open(writerObj);
for t_idx = 1:nt
    frame = M(t_idx);
    writeVideo(writerObj, M(t_idx));
end
close(writerObj)



  
