%% Parâmetros do problema
L = 1;                      % Comprimento do poço infinito
N = 1000;                   % Número de termos na série de Fourier
T = 100;                    % Tempo total de evolução
dt = 0.001;                 % Passo de tempo
hb = 1.05 * 10^-34;         % Constante de Planck reduzida
m = 9.109 * 10^-31;         % Massa da particula

%% Discretização espacial e temporal
dx = L/1000;                % Passo espacial
x = 0:dx:L;                 % Vetor de posições
t = 0:dt:T;                 % Vetor de tempos


%% Funções de Onda Inicial

% Distribuição gaussiana centrada no meio do poço)
mu = L/2;                                                       % Média da gaussiana
sigma = L/10;                                                   % Desvio padrão da gaussiana
gauss = exp(-(x-mu).^2/(2*sigma^2)) / sqrt(sigma*sqrt(2*pi));   % Função da Gaussiana

% Onda Degrau 
% x = linspace(0, L, N);                        % Vetor espacial
% dx = x(2) - x(1);                             % Passo espacial
% degrau = zeros(1, N);                         % Função de onda inicial zerada
% degrau(x < L/2) = sqrt(2/L);                  % Função da onda degrau

% Escolha de qual onda será utilizada
psi_0 = gauss;                                  % Aterar aqui para mudar qual onda será utilizada

% Normalização da função de onda
normFactor = sqrt(sum(abs(psi_0).^2)*dx);
psi_0 = psi_0 / normFactor;

%% Cálculo da evolução

% Cálculo dos coeficientes Cn
Cn = zeros(N, 1);
En = zeros(N, 1);
for n = 1:N
    En(n) = (n^2 * pi^2 * hb^2)/(2*m*L^2);      % Energia do n-ésimo estado
    Cn(n) = L * trapz(psi_0.*sin(n*pi*x/L),x);      % Produto interno com a função de onda inicial
end
    
% Evolução temporal
numSteps = ceil(T/dt); % Número de passos de tempo

ymax = 0;
for t = 0:N
    % Calculo da função de onda no determinado tempo
    psi = 0;
    for n = 1:N
        psi = psi + Cn(n)*exp(-1i*En(n)*t/hb).*sin(n*pi*x/L);
    end
        
    % Normalização da função de onda
    normFactor = sqrt(sum(abs(psi).^2)*dx);
    psi = psi / normFactor;

    % Define o valor maximo do eixo no primeiro loop
    if t == 0
        ymax = max(abs(psi).^2)*1.1;
    end
    
    % Visualização da curva de distribuição de probabilidade
    figure(1);
    plot(x, abs(psi).^2);
    xlabel('Posição');
    ylabel('Distribuição de Probabilidade');
    title(['Passo de Tempo: ', num2str(t)]);
    ylim([0, ymax]);
    drawnow;
end

