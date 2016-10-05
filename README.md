# TSL
Este programa aplica a Teoria Secular Linear para sistemas planetários.
Como trata-se de uma aproximação, consideramos pequenas inclinações e
baixas excentricidades.

Será usado o exemplo extraído de Caparpino, página 25, cujo arquivo pdf pode ser baixado [neste endereço](http://www.brera.mi.astro.it/~carpino/didattica/lagrange.pdf), para se verificar se o programa está sendo executado corretamente.

## Arquivos
* [TSL-Carpino.py](https://github.com/DeSouzaSR/TSL/blob/master/TSL-Carpino.py): Programa principal. Resolve a Teoria Secular Linear para um sistema planetário.
* [planetas-Carpino.txt](https://github.com/DeSouzaSR/TSL/blob/master/planetas-Carpino.txt): Dados de entrada. Cada linha do arquivo é um planeta.
* [Resultados-Carpino.txt](https://github.com/DeSouzaSR/TSL/blob/master/Resultados-Carpino.txt): Resultados da última execução do programa.

## Pré-requisitos
* Python 3.4.5
* Pandas 0.16.2
* Numpy 1.10.4
* Scipy 0.17.0
