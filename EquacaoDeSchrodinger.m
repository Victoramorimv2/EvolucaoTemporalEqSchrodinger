%% Par�metros do problema
L = 1;                      % Comprimento do po�o infinito
N = 1000;                   % N�mero de termos na s�rie de Fourier
T = 100;                    % Tempo total de evolu��o
dt = 0.001;                 % Passo de tempo
hb = 1.05 * 10^-34;         % Constante de Planck reduzida
m = 9.109 * 10^-31;         % Massa da particula

%% Discretiza��o espacial e temporal
dx = L/1000;                % Passo espacial
x = 0:dx:L;                 % Vetor de posi��es
t = 0:dt:T;                 % Vetor de tempos


%% Fun��es de Onda Inicial

% Distribui��o gaussiana centrada no meio do po�o)
mu = L/2;                                                       % M�dia da gaussiana
sigma = L/10;                                                   % Desvio padr�o da gaussiana
gauss = exp(-(x-mu).^2/(2*sigma^2)) / sqrt(sigma*sqrt(2*pi));   % Fun��o da Gaussiana

% Onda Degrau 
% x = linspace(0, L, N);                        % Vetor espacial
% dx = x(2) - x(1);                             % Passo espacial
% degrau = zeros(1, N);                         % Fun��o de onda inicial zerada
% degrau(x < L/2) = sqrt(2/L);                  % Fun��o da onda degrau

% Escolha de qual onda ser� utilizada
psi_0 = gauss;                                  % Aterar aqui para mudar qual onda ser� utilizada

% Normaliza��o da fun��o de onda
normFactor = sqrt(sum(abs(psi_0).^2)*dx);
psi_0 = psi_0 / normFactor;

%% C�lculo da evolu��o

% C�lculo dos coeficientes Cn
Cn = zeros(N, 1);
En = zeros(N, 1);
for n = 1:N
    En(n) = (n^2 * pi^2 * hb^2)/(2*m*L^2);      % Energia do n-�simo estado
    Cn(n) = trapz(psi_0.*sin(n*pi*x/L),x);      % Produto interno com a fun��o de onda inicial
end
    
% Evolu��o temporal
numSteps = ceil(T/dt); % N�mero de passos de tempo
psi = psi_0;

ymax = 0;
for t = 0:N
    % Calculo da fun��o de onda no determinado tempo
    for n = 1:N
        psi = psi + Cn(n)*exp(-1i*En(n)*t/hb).*sin(n*pi*x/L);
    end
        
    % Normaliza��o da fun��o de onda
    normFactor = sqrt(sum(abs(psi).^2)*dx);
    psi = psi / normFactor;

    % Define o valor maximo do eixo no primeiro loop
    if t == 0
        ymax = max(abs(psi).^2)*1.1;
    end
    
    % Visualiza��o da curva de distribui��o de probabilidade
    figure(1);
    plot(x, abs(psi).^2);
    xlabel('Posi��o');
    ylabel('Distribui��o de Probabilidade');
    title(['Passo de Tempo: ', num2str(t)]);
    ylim([0, ymax]);
    drawnow;
end

