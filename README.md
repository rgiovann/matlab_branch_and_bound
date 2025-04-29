# Algoritmo Branch and Bound para o Problema do Caixeiro Viajante (PCV)

Este projeto implementa uma solução para o **Problema de Roteirização** (também conhecido como **Problema do Caixeiro Viajante - TSP**) utilizando uma combinação de:
- **Algoritmo de Designação (Assignment Problem)**
- **Técnica de Branch and Bound**
- **Solver IBM CPLEX**, integrado via API ao MATLAB.

## Objetivo
Minimizar a distância total percorrida em uma rota que passa por todas as cidades apenas uma vez, retornando ao ponto de origem.

## Estrutura do Projeto
- `designacao_branch_and_bound.m`: Implementa o algoritmo principal de designação e branch and bound.
- `ChamaSolverCPLEX.m`: Função que chama o solver CPLEX para resolver as subproblemáticas de otimização.
- `branch_bound_roteirizacao.xlsx`: Exemplo de entrada de dados com as distâncias entre os pontos.

## Requisitos
- MATLAB instalado.
- Acesso ao **IBM CPLEX Optimizer** configurado para ser chamado via API no MATLAB.

## Como executar
1. Garanta que o CPLEX esteja instalado e integrado ao MATLAB.
2. Abra o MATLAB na pasta do projeto.
3. Execute o script `designacao_branch_and_bound.m`.
4. O programa lerá a matriz de distâncias (`branch_bound_roteirizacao.xlsx`) e buscará a rota ótima.

## Observações
- O método branch and bound é usado para reduzir o espaço de busca, descartando soluções que não podem melhorar o custo atual.
- O CPLEX é utilizado para resolver subproblemas de designação de maneira eficiente.
- O exemplo incluído é para estudos acadêmicos e pode ser expandido para casos maiores.

## Referências
- [Branch and Bound - Wikipedia](https://en.wikipedia.org/wiki/Branch_and_bound)
- [Problema do Caixeiro Viajante - Wikipedia](https://pt.wikipedia.org/wiki/Problema_do_caixeiro_viajante)