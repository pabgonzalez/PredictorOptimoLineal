%% Asegurarse de que "Archivo_2.mat" esté en la misma carpeta
fname='Archivo_2.mat';
data = load(fname);
x_size=data.ans(2);
N=128;

%% 1) Estimación de la autocorrelación total
varianz=var(data.x);
Rxx_np=zeros(1, N);
Rxx_p=zeros(1, N);
for k=1:N
    suma=0;
    for i=1:(x_size-k)
        suma = suma + data.x(i)*data.x(i+k-1);
    end
    Rxx_np(k)=(1/(x_size-k+1))*suma/varianz;
    Rxx_p(k)=(1/x_size)*suma/varianz;
end

plot(1:128, Rxx_np, ' o ','MarkerSize',2)
xlim([0, 130])
xlabel('k')
ylabel('r_{XX}(k)')

hold on
plot(1:128, Rxx_p, ' o ','MarkerSize',3)

legend('No polarizado', 'Polarizado')

%% 2) Estimación de autocorrelación parcial
%Armo matriz de la eq de Yule-Walker

revRxx=fliplr(Rxx_np);
supreme_fila=[revRxx, 1, Rxx_np];
phis=zeros(1, N);

phis(1)=Rxx_np(1);
for orden=2:N
    matA=supreme_fila((N+1):(N+1+(orden-1)));     %quiero estimar Phi(orden, orden)
    for n=1:(orden-1)
        matA=[matA;supreme_fila((N+1-n):(N+1+orden-1-n))];
    end
    invA=inv(matA);
    pp=invA*transpose(Rxx_np(1:orden));
    phis(orden)= pp(orden);
end

figure
plot(1:N, phis, ' . ','MarkerSize',5)
xlim([0, 130])
xlabel('k')
ylabel('\phi_{k,k}')

%% 3) Parámetros del modelo
disp(['-----------------'])
disp(['Primeros 10 valores de rxx:'])
disp([num2str(Rxx_np(1:10))])
disp(['Se consideran 0 a partir de rxx(5), por lo que queda que es modelo MA(4).'])
disp([' '])

Th_calc=zeros(1, 4);
Th_calc(4)=Rxx_np(4);
Th_calc(3)= Rxx_np(3)/(Rxx_np(4)+1);
Th_calc(2)= (Rxx_np(2)-Th_calc(3)*Th_calc(4))/(Rxx_np(3)+1);
Th_calc(1)= (Rxx_np(1)-Th_calc(2)*Th_calc(3)-Th_calc(3)*Th_calc(4))/(Rxx_np(2)+1);

disp(['Si la entrada es una secuencia de ruido blanco y Gaussiano con varianza unitaria:'])
disp(['X(n) = (', num2str(Th_calc(1)), ')*e(n-1) + (', num2str(Th_calc(2)),')*e(n-2) + (', num2str(Th_calc(3)),')*e(n-3) + (', num2str(Th_calc(4)),')*e(n-4) + e(n)'])
disp([num2str(Th_calc)])

disp(['-----------------'])

%% 4) Cálculo analítico
varN=varianz/(1+sum(Th_calc.^2));

Rxx_teo=zeros(1, 128);
for k=1:3
    m=0;
    for j=(k+1):4
        m = m + Th_calc(j)*Th_calc(j-k);
    end
    Rxx_teo(k)=varN*(Th_calc(k)+m);
end
Rxx_teo(4)=varN*Th_calc(4);

figure
plot(1:128, Rxx_teo)
hold on
plot(1:128, Rxx_np)
xlabel('k')
ylabel('r_{XX}(k)')
legend('r_{XX} analítico', 'r_{XX} estimado')
%% 5) Densidad espectral de potencia

%Promediación de periodogramas
Xf=fftshift(fft(data.x));
Snn_p=(abs(Xf).^2)/x_size;

%DFT de Rxx
Snn=fftshift(fft(Rxx_np));

figure
plot((-(length(Snn_p)-1)/2:(length(Snn_p)-1)/2)*100/length(Snn_p), abs(Snn_p))

hold on
plot((-(length(Snn)-1)/2:(length(Snn)-1)/2)*100/length(Snn), abs(Snn))
xlabel('f')
ylabel('S_{NN}(f)')

legend('Promediación de periodogramas', 'DFT de Rxx')



