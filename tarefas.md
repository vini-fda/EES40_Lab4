# Tarefas

## Parte 1

Projeto DLQR da lei de realimentação ótima de estado discreto no tempo. O modelo contínuo SISO LIT tem função de transferência:

$$ G(s) = \frac{Y(s)}{\hat{U}(s)} = \frac{1000}{(s+10)(s+0.5+10j)(s+0.5-10j)} $$

### Tarefa 1

Desconsiderando-se os ruídos, com o modelo determinístico projetar a realimentação de estado discreto 

$$ u[k] = -K_\text{dlqr} \cdot x[k] $$

Verificar se a resposta amostrada à referência degrau da malha fechada do modelo determinístico discreto com saída $$[y \quad y']^T$$ com realimentação de estado ótima DLQR atende aos requisitos.

## Parte 2

Filtro de Kalman discreto no tempo e controle DLQG

### Tarefa 2

Sintonize o filtro de Kalman discreto apresentado no _script_ de projeto, simule o controlador DLQG; use a função de autocorrelação normalizada para avaliar a hipótese de brancura da sequência de inovação $$ y'[k] - C_{k} \hat{x}_\text{aug} [k|k-1] $$ segundo os resultados da simulação. Em seguida, use o _script_ que controla a placa AD/DA da NI, `EES40_2023_PlacaNlLab4controller_DLQGcomKF.m`, para implementar o controlador DLQG na bancada. A frequência máxima de amostragem do D/A é 150 Hz. Verifique os resultados e analise se a hipótese é válida.

## Tarefa 3

Projetar filtro de Kalman de ordem três que estima somente o estado da planta com a medida da saída da planta \( y \) e o empregar no controlador DLQG. Vide o _script_ de projeto `Lab4_EES40_2023_semana14_DLQGcomKFredux_alunado.m`. Implementar na bancada com o _script_ `EES40_2023_PlacaNI_Lab4controller_DLQGcomKFredux.m` e sintonizar. A hipótese de brancura da sequência de inovação se sustenta na operação em bancada? Explique.
