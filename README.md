# ⚡ Ejemplo Básico de Detección de Eventos en Cargas Eléctricas Residenciales (NILM)

Este repositorio contiene la implementación del algoritmo **NILM** (Non-Intrusive Load Monitoring) para la detección de eventos de encendido y apagado de un refrigerador (heladera), a partir de una única señal de consumo eléctrico. El método utiliza la **Transformada Wavelet Estacionaria (SWT)** y un **Criterio Dual** para identificar de forma robusta los cambios transitorios de potencia.

## 📝 Metodología

El proyecto se basa en el **pipeline** metodológico descrito en el artículo "[Detección de Eventos en Cargas Eléctricas Residenciales Mediante Transformada Wavelet Estacionaria y Criterio Dual](AEAndradaTivani_DetecciondeEventosenCargasElectricasResidencialesSWTyCriterio Dual.pdf)". La metodología consta de los siguientes pasos:

1.  **Preprocesamiento de datos**: Carga y limpieza de una serie temporal de consumo de potencia (VA).
2.  **Descomposición Wavelet**: Aplicación de la SWT para descomponer la señal en diferentes escalas de tiempo-frecuencia, aislando los eventos transitorios del consumo de fondo y el ruido.
3.  **Análisis de Energía**: Cálculo de la energía de los coeficientes wavelet en los niveles de detalle relevantes.
4.  **Detección de Picos**: Uso de un umbral adaptativo y robusto, basado en la Mediana de la Desviación Absoluta (**MAD**), para detectar picos de energía significativos.
5.  **Criterio Dual**: Confirmación de los eventos detectados mediante un segundo criterio sobre la **magnitud del cambio de potencia ($\Delta P$)** para descartar falsos positivos.
6.  **Fusión de Eventos**: Unificación de las detecciones para generar una tabla final con los tiempos de inicio y fin de cada ciclo de evento.

Este enfoque evita el uso de hardware especializado y **wavelets a medida**, lo que lo hace más versátil y escalable para diferentes aplicaciones.

## 📂 Estructura del Repositorio

-   `NILM_Residential_detectionWavelet.m`: Script principal de MATLAB que implementa la metodología descrita.
-   `ESP32Power - PowerLogs(14,15 del 09).csv`: Archivo de ejemplo con datos de consumo eléctrico residencial (Vrms, Irms) utilizados para las pruebas.
-   `AEAndradaTivani_DetecciondeEventosenCargasElectricasResidencialesSWTyCriterio Dual.pdf`: Tesis o artículo científico que fundamenta el desarrollo del proyecto.

## ⚙️ Requisitos

-   MATLAB (se recomienda la versión R2020b o superior).
-   Se puede requerir la Toolbox de Wavelet de MATLAB.

## 🚀 Uso

1.  Clona este repositorio en tu máquina local.
2.  Abre el archivo `NILM_Residential_detectionWavelet.m` en MATLAB.
3.  Asegúrate de que los archivos de datos `.csv` se encuentren en la misma carpeta o actualiza la ruta de carga en el script.
4.  Ejecuta el script. El programa procesará los datos, detectará los eventos y generará visualizaciones de los resultados.

---

## 🤝 Agradecimientos

Agradecemos especialmente al **Dr. Ing. Juan Astrada** por su invaluable colaboración en la obtención de los datos de consumo eléctrico utilizados en este proyecto.

## 📜 Licencia

Este proyecto está bajo la Licencia MIT. Consulta el archivo `LICENSE` para más detalles.

## 📧 Contacto

Para cualquier pregunta o sugerencia, puedes contactar a los autores:

-   **Esp. Ing. Astri Edith Andrada Tivani**: astriandrada@ing.unrc.edu.ar
-   **Dr. Ing. Juan Astrada**: jastrada@ing.unrc.edu.ar
