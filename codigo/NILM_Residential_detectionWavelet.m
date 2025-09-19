%% ALGORITMO DE DETECCIÓN DE EVENTOS EN CARGAS ELÉCTRICAS (NILM) USANDO DWT
%
% OBJETIVO:
% Identificar automáticamente los ciclos de encendido/apagado de un
% dispsitivo residencial (ej. una heladera) a partir de la señal de consumo 
% eléctrico total de una vivienda.
%
% METODOLOGÍA:
% 1.  Carga y preprocesamiento de una serie temporal de consumo de potencia.
% 2.  Aplicación de la Transformada Wavelet Estacionaria (SWT) para 
%     descomponer la señal en diferentes escalas de tiempo/frecuencia. 
%     Esto permite separar los cambios bruscos (eventos) del ruido y del 
%     consumo de fondo.
% 3.  Calculo de la energía de los coeficientes wavelet en los niveles de
%     interés, donde se espera encontrar las "firmas" de los eventos.
% 4.  Definición de un umbral adaptativo y robusto, basado en la Mediana de 
%     la Desviación Absoluta (MAD), para detectar picos de energía.
% 5.  Aplicación de un segundo criterio sobre la magnitud del cambio de 
%     potencia (ΔP) para confirmar que el evento es significativo.
% 6.  Limpiieza y fusión de las detecciones para generar una tabla de 
%     eventos final, con su hora de inicio, fin y duración.
%
% AUTOR: Astri Edith Andrada Tivani
% FECHA: 17 de septiembre de 2025
% --------------------------------------------------------------------------
%% Comandos Iniciales de Entorno
% --------------------------------------------------------------------------
clear;          % Limpia el workspace de variables
close all;      % Cierra todas las figuras
clc;            % Limpia la ventana de comandos
fprintf('Iniciando algoritmo de detección NILM basado en DWT...\n');
% --------------------------------------------------------------------------
%% SECCIÓN 1: PARÁMETROS DE CONFIGURACIÓN
% --------------------------------------------------------------------------
fprintf('1) Configurando parámetros...\n');
% --- Parámetros de Entrada de Datos ---
% Archivo con los datos de potencia.
filename = 'ESP32Power - PowerLogs(14,15 del 09).csv'; 
% Factor de corrección para el sensor de corriente.
I_gain = 10;               
% --- Parámetros del Análisis Wavelet ---
% La 'Daubechies 4' es una wavelet comúnmente usada en análisis de señales
% de potencia por su buen balance entre suavidad y capacidad de localización.
wavelet_name = 'db4';       
% Nivel de descomposición. Un nivel 5 permite analizar la señal en
% escalas de tiempo que van desde muy rápidas (ruido) a más lentas,
% donde se espera encontrar la firma de los dispositivos residenciales.
decomp_level = 5;
% Se utiliza la Transformada Wavelet Estacionaria (SWT) en lugar de la DWT
% estándar. La SWT es invariante a la traslación, lo que significa que 
% detecta un evento sin importar en qué punto exacto de la señal comience.
use_swt = true;     
% --- Parámetros de Detección de Eventos ---
% Factor de sensibilidad del umbral. Un valor más alto hace la detección más
% estricta (menos falsos positivos). Este valor se determinó empíricamente.
mad_k = 4.8; 
% Umbral mínimo de cambio de potencia para considerar un evento válido.
% Evita que picos de energía wavelet causados por ruido sean detectados.
deltaP_threshold_VA = 80;   % [VA]
% --- Parámetros de Visualización ---
plot_figures = true;        % Controla si se generan los gráficos de resultados.
% --------------------------------------------------------------------------
%% SECCIÓN 2: CARGA Y PREPARACIÓN DE LOS DATOS
% --------------------------------------------------------------------------
fprintf('2) Cargando y preparando los datos desde ''%s''...\n', filename);
% 2.1) Lectura de datos y parseo de timestamps
opts = detectImportOptions(filename);
T = readtable(filename, opts);
% Se busca un formato de fecha y hora adecuado, probando varios formatos comunes
% para asegurar la compatibilidad con diferentes archivos de origen.
if any(ismember(T.Properties.VariableNames, {'Date','Time'}))
    timeStrings = string(T.Date) + " " + string(T.Time);
    timeStrings = replace(timeStrings, ["a. m.", "p. m."], ["AM", "PM"]); % Estandarizar AM/PM
    
    % --- BLOQUE DE LECTURA DE FECHA ROBUSTO ---
    try
        % INTENTO 1: Formato con mes como texto y AM/PM (ej. 17-Sep-2025 09:21:10 AM)
        timestamps_full = datetime(timeStrings, 'InputFormat', 'dd-MMM-yyyy hh:mm:ss a');
    catch ME1
        fprintf('Formato ''dd-MMM-yyyy...'' no funcionó. Intentando otro...\n');
        try
            % INTENTO 2: Formato numérico estándar (ej. 17/09/2025 21:21:10)
            timestamps_full = datetime(timeStrings, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
        catch ME2
            fprintf('Formato ''dd/MM/yyyy...'' tampoco funcionó. Intentando detección automática...\n');
            try
                % INTENTO 3: Dejar que MATLAB reconozca el formato.
                timestamps_full = datetime(timeStrings);
            catch ME3
                % Si todo falla, muestra un error claro y un ejemplo del texto problemático.
                disp('ERROR: No se pudo convertir el texto a fecha y hora con ninguno de los métodos.');
                disp('Por favor, revisa el formato de la primera línea de tu archivo CSV:');
                disp(timeStrings(1));
                rethrow(ME3); % Re-lanza el último error para detener la ejecución.
            end
        end
    end
    % Es crucial ordenar los datos cronológicamente por si el registro no fue secuencial.
    [timestamps_full, sort_idx] = sort(timestamps_full);
    T = T(sort_idx,:);
else
    timestamps_full = T.Timestamp;
end
% 2.2) Limpieza y cálculo de la Potencia Aparente
% Se eliminan lecturas inválidas o con errores (ej. Vrms o Irms = 0).
idx_valid = T.Vrms > 0 & T.Irms > 0;
T_full = T(idx_valid,:);
timestamps_full = timestamps_full(idx_valid);
% Se aplica el factor de ganancia para obtener la corriente real.
Irms_full = T_full.Irms / I_gain;
Vrms_full = T_full.Vrms;
% Se calcula la Potencia Aparente (S), que es el producto de Vrms e Irms.
% Esta es la señal principal que se analizará.
PowerVA_full = Vrms_full .* Irms_full;
% --------------------------------------------------------------------------
%% SECCIÓN 3: VISUALIZACIÓN INICIAL DE LA SEÑAL COMPLETA
% --------------------------------------------------------------------------
% Este gráfico ofrece una visión general de todo el conjunto de datos y
% resalta visualmente la ventana de tiempo que será analizada en detalle.
if plot_figures
    fprintf('3) Generando gráfico de la señal original completa...\n');
    figure('Name', 'Visión General - Señal Original Completa');
    hold on;
    plot(timestamps_full, PowerVA_full, 'k-', 'DisplayName', 'Consumo Total (Original)');
    
    % Se define la ventana de análisis nocturna (22:30 a 06:00).
    % Este horario se elige para minimizar la interferencia de otros aparatos.
    h_full = hour(timestamps_full);
    m_full = minute(timestamps_full);
    idx_noche_full = ((h_full == 22 & m_full >= 30) | h_full > 22 | h_full < 6);
    
    % Se dibuja un rectángulo sombreado para marcar la región de interés.
    start_idx = find(idx_noche_full, 1, 'first');
    end_idx = find(idx_noche_full, 1, 'last');
    if ~isempty(start_idx) && ~isempty(end_idx)
        patch_x = [timestamps_full(start_idx), timestamps_full(end_idx), timestamps_full(end_idx), timestamps_full(start_idx)];
        ylim_vals = ylim;
        patch_y = [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)];
        patch(patch_x, patch_y, 'red', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Ventana de Análisis (22:30-06:00h)');
    end
    title('Señal Original y Ventana de Análisis');
    ylabel('Potencia Aparente (VA)'); xlabel('Timestamp'); legend('show'); grid on; hold off;
end
% --------------------------------------------------------------------------
%% SECCIÓN 4: PREPROCESAMIENTO DE LA SEÑAL DE TRABAJO
% --------------------------------------------------------------------------
fprintf('4) Preprocesando la señal para el análisis...\n');
% 4.1) Estimación de la frecuencia de muestreo (fs)
% Se usa la mediana para robustez ante muestras irregulares
dt = median(seconds(diff(timestamps_full))); 
% Fallback por si hay un solo dato
if isnan(dt) || dt <= 0, dt = 1; end 
fs = 1/dt;
% 4.2) Recorte a la ventana de análisis (22:30 - 06:00)
h = hour(timestamps_full); 
m = minute(timestamps_full);
idx_noche = ((h == 22 & m >= 30) | h > 22 | h < 6);
timestamps = timestamps_full(idx_noche);
PowerVA = PowerVA_full(idx_noche);
% 4.3) Cálculo de la línea base (Baseline) para eliminar la tendencia
% Se usa una media móvil basada en la mediana (movmedian) sobre una ventana larga
% (2 horas) para estimar el consumo de fondo. La mediana es robusta a los
% picos de consumo de los eventos que queremos detectar.
% 7200 segundos = 2 horas
baseline = movmedian(PowerVA, round(fs * 7200)); 
% 4.4) Extracción de la señal de detalle (Detrending)
% Al restar la línea base, la señal resultante (PowerVA_detr) se centra en cero.
% Ahora, cualquier valor > 0 representa un consumo por encima de lo normal,
% es decir, un potencial evento.
PowerVA_detr = PowerVA - baseline;
% 4.5) Imputación de datos faltantes (si los hubiera)
% La transformada wavelet es sensible a discontinuidades. Se rellenan posibles
% valores NaN con una mediana móvil corta para no introducir artefactos.
PowerVA_detr = fillmissing(PowerVA_detr, 'movmedian', round(fs * 60));
% La señal final 'PowerVA_clean' es la que se emplea para el análisis wavelet.
PowerVA_clean = PowerVA_detr;
% --------------------------------------------------------------------------
%% SECCIÓN 5: ANÁLISIS MULTI-RESOLUCIÓN CON TRANSFORMADA WAVELET
% --------------------------------------------------------------------------
fprintf('5) Aplicando la Transformada Wavelet Estacionaria (SWT)...\n');
% 5.1) Padding para compatibilidad con SWT
% La SWT requiere que la longitud de la señal sea divisible por 2^L, donde L
% es el nivel de descomposición. Se añade un relleno simétrico (padding)
% al final de la señal para cumplir este requisito.
n = length(PowerVA_clean);
fprintf('Longitud original de la señal: n=%d\n', n);
if use_swt
    L = decomp_level;
    required_divisor = 2^L;
    new_len = ceil(n / required_divisor) * required_divisor;
    fprintf('Rellenando señal de %d a %d muestras para compatibilidad con SWT.\n', n, new_len);
    signal_padded = wextend('1D', 'sym', PowerVA_clean, new_len - n, 'r');
    
    % 5.2) Aplicación de la SWT
    swtC = swt(signal_padded, L, wavelet_name);
    
    % 5.3) Extracción de Coeficientes de Detalle
    % Se extraen los coeficientes de detalle D1..D5. Cada nivel 'D_j' captura
    % eventos en una escala de tiempo diferente:
    % - D1: Cambios muy rápidos, usualmente ruido.
    % - D2-D4: Transitorios de corta y media duración, ideal para ON/OFF de aparatos.
    % - D5: Variaciones más lentas en la señal.
    D = swtC(1:L, 1:n); % Se descarta el padding
else
    % (Código alternativo para DWT estándar)
    [C,L_dwt] = wavedec(PowerVA_clean, decomp_level, wavelet_name);
    D = zeros(decomp_level, n);
    for lev = 1:decomp_level
        D(lev,:) = wrcoef('d', C, L_dwt, wavelet_name, lev);
    end
end
fprintf('Descomposición completada.\n');
% --------------------------------------------------------------------------
%% SECCIÓN 6: DETECCIÓN Y PROCESAMIENTO DE EVENTOS
% --------------------------------------------------------------------------
fprintf('6) Detectando eventos a partir de los coeficientes wavelet...\n');
% 6.1) Criterio 1: Umbral de Energía Wavelet
% Se seleccionan los niveles de detalle que capturan mejor los eventos de interés.
% Empíricamente, los niveles 3 y 4 son efectivos para dispositivos residenciales comunes.
levels_of_interest = intersect([3 4], 1:decomp_level);
% Se calcula la energía instantánea como el cuadrado de los coeficientes (D^2).
energy_per_level = D.^2;
% Se suaviza la señal de energía con una media móvil para obtener una "envolvente".
wlen = max(1, round(12 * fs)); % Ventana de suavizado
energy_env = movmean(energy_per_level', wlen)';
% Se suma la energía de los niveles de interés para obtener una única señal de energía.
tot_energy = sum(energy_env(levels_of_interest, :), 1)';
% Se calcula un umbral estadístico robusto. Se usa la Mediana (med) y la
% Desviación Absoluta Mediana (MAD) porque no son afectadas por los propios
% eventos (que son 'outliers' en la señal de energía).
med = median(tot_energy);
madv = mad(tot_energy, 1); % El segundo argumento '1' es para compatibilidad con versiones antiguas
threshold = med + mad_k * madv;
fprintf('Umbral de energía = %.3f (med=%.3f, mad=%.3f, k=%.1f)\n', threshold, med, madv, mad_k);
% Se crea una máscara booleana con los puntos que superan el umbral.
mask_energy = tot_energy > threshold;
% 6.2) Criterio 2: Umbral de Cambio de Potencia (ΔP)
% Este es un criterio de validación. Un evento real no solo debe generar un
% pico de energía wavelet, sino también un cambio de potencia significativo.
% Se calcula ΔP restando una línea base local de corto plazo.
win_base_local = max(1, round(30 * fs));
local_baseline = movmedian(PowerVA, win_base_local);
deltaP = PowerVA - local_baseline;
% Se crea una segunda máscara para los puntos que superan el umbral de ΔP.
mask_deltaP = deltaP > deltaP_threshold_VA;
% 6.3) Combinación de Criterios y Limpieza Morfológica
% Un evento válido debe cumplir AMBOS criterios (AND lógico).
mask_combined = mask_energy & mask_deltaP;
% Se realiza una operación de "cierre morfológico" para unir segmentos de
% detección que estén muy cercanos, formando un único evento continuo.
se_len = max(1, round(5 * fs));
mask_closed = conv(double(mask_combined), ones(1, se_len), 'same') > 0;
% Se eliminan detecciones muy cortas que probablemente sean ruido.
min_run = round(3 * fs);
cc_tmp = bwconncomp(mask_closed);
mask_clean = false(size(mask_closed));
for ii = 1:cc_tmp.NumObjects
    idxs = cc_tmp.PixelIdxList{ii};
    if numel(idxs) >= min_run
        mask_clean(idxs) = true;
    end
end
% 6.4) Estructuración y Fusión de Eventos
% Se convierte la máscara booleana final en una tabla de eventos.
min_event_duration_s = 8; % Duración mínima en segundos para ser un evento válido
min_samples = round(min_event_duration_s * fs);
cc = bwconncomp(mask_clean);
event_table = table();
ev_idx = 0;
for k = 1:cc.NumObjects
    idxs = cc.PixelIdxList{k};
    if numel(idxs) >= min_samples
        ev_idx = ev_idx + 1;
        event_table(ev_idx,:) = {ev_idx, timestamps(min(idxs)), timestamps(max(idxs)), numel(idxs)/fs, max(tot_energy(idxs))};
    end
end
% Se definen los nombres de las columnas de la tabla.
if isempty(event_table)
    event_table = cell2table(cell(0,5), 'VariableNames', {'EventID','StartTime','EndTime','Duration_s','PeakEnergy'});
else
    event_table.Properties.VariableNames = {'EventID','StartTime','EndTime','Duration_s','PeakEnergy'};
end
% A veces, un solo ciclo de un aparato puede ser detectado como dos eventos
% muy seguidos. Se fusionan eventos que tengan un espacio menor a `max_gap_s`.
max_gap_s = 90; % [segundos]
if height(event_table) > 1
    newtab = []; i = 1;
    while i <= height(event_table)
        st = event_table.StartTime(i);
        ed = event_table.EndTime(i);
        j = i + 1;
        while j <= height(event_table) && seconds(event_table.StartTime(j) - ed) <= max_gap_s
            ed = event_table.EndTime(j);
            j = j+1;
        end
        idxs = find(timestamps >= st & timestamps <= ed);
        if isempty(idxs)
            pk = NaN;
            dur = seconds(ed - st);
        else
            pk = max(tot_energy(idxs));
            dur = numel(idxs) / fs;
        end
        newtab = [newtab; {size(newtab,1)+1, st, ed, dur, pk}];
        i = j;
    end
    event_table = cell2table(newtab, 'VariableNames', {'EventID','StartTime','EndTime','Duration_s','PeakEnergy'});
end
% --------------------------------------------------------------------------
%% SECCIÓN 7: VISUALIZACIÓN DE RESULTADOS Y SALIDA
% --------------------------------------------------------------------------
fprintf('7) Generando gráficos de resultados y mostrando tabla final...\n');
if plot_figures
    % Llamada a las funciones encargadas de generar los gráficos.
    PlotResults(timestamps, PowerVA, baseline, tot_energy, threshold, D, energy_env, event_table);
    if height(event_table) > 0
        PlotDetailedAnalysis(timestamps, PowerVA, PowerVA_clean, wavelet_name, event_table, tot_energy, threshold, deltaP, deltaP_threshold_VA);
        PlotEnergyAnalysis(timestamps, PowerVA, energy_env, event_table);
    end
end
% Se muestra el resultado final en la consola.
disp('Ciclos de heladera detectados:');
disp(event_table);
fprintf('Proceso finalizado.\n');
% --------------------------------------------------------------------------
%% SECCIÓN 8: FUNCIONES AUXILIARES DE PLOTEO
% --------------------------------------------------------------------------
%% PlotResults
% Genera una vista multinivel del proceso de detección completo.
% Muestra la señal original, la energía combinada con su umbral, y la
% descomposición wavelet nivel por nivel, resaltando los eventos detectados.
function PlotResults(timestamps, PowerVA, baseline, tot_energy, threshold, D, energy_env, event_table)
    decomp_level = size(D,1);
    figure('Name','Resultados de Detección Final','Units','normalized','Position',[0.03 0.03 0.94 0.85]);
    ax = gobjects(decomp_level + 2, 1);
    
    % Subplot 1: Señal de potencia original y su baseline
    ax(1) = subplot(decomp_level + 2, 1, 1);
    plot(timestamps, PowerVA, 'b'); hold on; plot(timestamps, baseline, 'k--');
    ylabel('Potencia (VA)'); title('Señal de Trabajo (22:30 - 06:00)'); legend('Potencia','Baseline'); grid on;
    
    % Subplot 2: Energía combinada y umbral de detección
    ax(2) = subplot(decomp_level + 2, 1, 2);
    plot(timestamps, tot_energy, 'k'); hold on; yline(threshold,'r--','Umbral');
    ylabel('Energía Wavelet'); title('Energía Combinada (D3+D4) y Umbral'); grid on;
    
    % Subplots 3 a N: Coeficientes de detalle y su energía por nivel
    for lev=1:decomp_level
        ax(2+lev) = subplot(decomp_level + 2, 1, 2+lev);
        plot(timestamps, D(lev,:), 'r'); hold on;
        yyaxis right;
        plot(timestamps, energy_env(lev,:), 'k--');
        ylabel('Energía'); gca.YAxis(1).Color = 'r';
        title(sprintf('Coeficientes de Detalle D%d y su Envolvente de Energía', lev));
    end
    
    linkaxes(ax,'x'); xlim([timestamps(1) timestamps(end)]);
    
    % Resaltar los eventos detectados en todos los subplots
    % >>> MODIFICACIÓN: remarcar áreas verdes <<<
    for i=1:height(event_table)
        for k=1:numel(ax)
            axes(ax(k)); hold on;
            patch_x = [event_table.StartTime(i), event_table.EndTime(i), event_table.EndTime(i), event_table.StartTime(i)];
            ylim_i = ylim; 
            patch_y = [ylim_i(1), ylim_i(1), ylim_i(2), ylim_i(2)];
            patch(patch_x, patch_y, 'green','FaceAlpha',0.4,'EdgeColor','green','LineWidth',1.2);
        end
    end
    annotation('textbox',[0.8 0.1 0.1 0.1], 'String', sprintf('Ciclos Detectados: %d', height(event_table)),'FitBoxToText','on','BackgroundColor','white');
end
%% PlotDetailedAnalysis
% Ofrece una vista detallada de los pasos clave del preprocesamiento y
% de los dos criterios de detección, haciendo zoom en el primer evento.
function PlotDetailedAnalysis(timestamps, PowerVA, PowerVA_clean, wavelet_name, event_table, tot_energy, threshold, deltaP, deltaP_threshold_VA)
    figure('Name', 'Análisis Detallado de Criterios', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    
    % Gráfico 1: Comparación de la señal antes y después del preprocesamiento
    subplot(2, 2, 1);
    plot(timestamps, PowerVA, 'b-', 'DisplayName', '1. Señal Recortada'); hold on;
    plot(timestamps, PowerVA_clean, 'r-', 'DisplayName', '2. Señal de Trabajo (sin tendencia)');
    title('Paso 1: Preprocesamiento de la Señal'); legend; grid on; ylabel('Potencia (VA)');
    
    % Gráfico 2: Visualización de la wavelet madre utilizada
    subplot(2, 2, 2);
    try
        [~, psi, ~] = wavefun(wavelet_name, 10); plot(psi, 'LineWidth', 2); title(['Paso 2: Forma de la Wavelet Madre (' wavelet_name ')']); grid on;
    catch
        title('Error al graficar la wavelet');
    end
    
    % Se hace zoom en el primer evento detectado para ver los detalles
    primer_evento = event_table(1,:);
    margen = seconds(60); 
    idx_zoom = find(timestamps >= (primer_evento.StartTime - margen) & timestamps <= (primer_evento.EndTime + margen));
    
    % Gráfico 3: Aplicación del primer criterio (Energía Wavelet)
    subplot(2, 2, 3);
    yyaxis left; plot(timestamps(idx_zoom), PowerVA(idx_zoom),'-','LineWidth',1.5,'DisplayName','Potencia'); ylabel('Potencia (VA)');
    yyaxis right; plot(timestamps(idx_zoom), tot_energy(idx_zoom),'-','LineWidth',1.5,'DisplayName','Energía'); hold on; yline(threshold, 'r--', 'DisplayName', 'Umbral Energía');
    title('Paso 3: Criterio #1 (Energía Wavelet)'); legend; grid on;
    
    % Gráfico 4: Aplicación del segundo criterio (Cambio de Potencia ΔP)
    subplot(2, 2, 4);
    yyaxis left; plot(timestamps(idx_zoom), PowerVA(idx_zoom),'-','LineWidth',1.5,'DisplayName','Potencia'); ylabel('Potencia (VA)');
    yyaxis right; plot(timestamps(idx_zoom), deltaP(idx_zoom),'-','LineWidth',1.5,'DisplayName','\DeltaP'); hold on; yline(deltaP_threshold_VA, '--', 'DisplayName', 'Umbral \DeltaP');
    title('Paso 4: Criterio #2 (Cambio de Potencia)'); legend; grid on;
end
%% PlotEnergyAnalysis
% Genera un gráfico similar a los que se ven en papers académicos, mostrando
% la energía de cada nivel de detalle superpuesta a la señal de potencia.
function PlotEnergyAnalysis(timestamps, PowerVA, energy_env, event_table)
    decomp_level = size(energy_env, 1);
    figure('Name', 'Análisis de Energía por Nivel', 'Units', 'normalized', 'Position', [0.2 0.05 0.6 0.85]);
    ax = gobjects(decomp_level, 1);
    
    for lev = 1:decomp_level
        ax(lev) = subplot(decomp_level, 1, lev);
        hold on;
        % Se dibuja la potencia de fondo como referencia visual.
        plot(timestamps, PowerVA, 'm--', 'DisplayName', 'Potencia (VA)');
       
        % Se dibuja la envolvente de energía de este nivel.
        plot(timestamps, energy_env(lev, :), 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('Energía D%d', lev));
        
        ylabel('Energía');
        title(sprintf('Energía en Nivel de Detalle D%d', lev));
        grid on;
        legend('show', 'Location', 'northeast');
    end
    
    % Se resaltan los eventos finales detectados en todos los niveles.
    % >>> MODIFICACIÓN: remarcar áreas verdes <<<
    for i=1:height(event_table)
        for k=1:numel(ax)
            axes(ax(k));
            hold on;
            patch_x = [event_table.StartTime(i), event_table.EndTime(i), event_table.EndTime(i), event_table.StartTime(i)];
            ylim_i = ylim;
            patch_y = [ylim_i(1), ylim_i(1), ylim_i(2), ylim_i(2)];
            patch(patch_x, patch_y, 'green','FaceAlpha',0.4,'EdgeColor','green','LineWidth',1.2);
        end
    end
    linkaxes(ax, 'x');
    xlim([timestamps(1) timestamps(end)]);
end
