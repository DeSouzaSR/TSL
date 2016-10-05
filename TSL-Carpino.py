# coding: utf-8

# # Teoria Secular Linear
# 
# O objetivo deste programa é aplicar a Teoria Secular Linear a um exemplo
# extraído de Carpino Teoria lineare delle perturbazioni secolari
#  (o texto está em italiano). O arquivo pode ser baixado
# neste site: http://www.brera.mi.astro.it/~carpino/didattica/lagrange.pdf. 
# 
# ## Definições
# * a: Semieixo maior
# * e: Excentricidade
# * I: Inclinação
# * capom: Longitude do nodo ascendente
# * omega: Argumento do periélio
# * M: Anomalia média
# * varpi = capom + omega: Longitude do periélio
# * lambda = M + varpi: Longitude média 

#
# ## Procedimento para obter a solução das equações
# * Calcular alpha_{ij} para todos os corpos envolvidos
# * Calcular os coeficientes de Laplace
# * Calcular os coeficientes das matrizes A e B
# * Calcular os autovetores e autovalores das matrizes simétricas A e B
# * Calcular os autovetores G_{ij} e S_{ij}
# * Determinar as condições iniciais Gamma_{i}, eta_{i}, Sigma, nu
#   a partir dos elementos médios


# %% Importando módulos
import numpy as np
import pandas as pd
from scipy.integrate import quad

# %% Definições de funções
# Função alpha
def al(a1, a2):
    """ Função que calcula a razão entre o mínimo e o máximo de dois semieixos
        maiores
        
        Input: a1, a2: dois semieixo maiores
        
        Output: min(a1, a2) / max(a1, a2)
    """
    return min(a1, a2) / max(a1, a2)

# %% Coeficientes de Laplace $b_{\psi}^{k}(\alpha_{i,j})$ 
def coeff_laplace(K, al):
    """ Função que calcula os coeficientes de Laplace"""
    f = lambda x: (1 - 2 * al * np.cos(x) + al**2)**(-3/2) * np.cos(K * x)
    coeff = quad(f, 0, 2*np.pi)
    return (1.0 / np.pi) * coeff[0]
    
# %% Definição das funções parciais para o cálculo das matrizes
def parte1(i, G, n, a):
    """ Função que calcula a equação G/(4 * n[i] * a[i]**2)
        como parte da determinação da matrizes A e B    
    """
    return G/(4 * n[i] * a[i]**2)

def parte2(i, j, m, al, a):
    """ Função que calcula a equação (m[j] * al(a[i], a[j])) / max(a[i], a[j])
        como parte da determinação da matrizes A e B    
    """
    return (m[j] * al(a[i], a[j])) / max(a[i], a[j])

def parte3(K, i, j, coeff_laplace, a):
    """ Função que calcula a equação coeff_laplace(K, al(a[i], a[j]))
        como parte da determinação da matrizes A e B    
    """
    return coeff_laplace(K, al(a[i], a[j]))

# %% Função que adiciona título com sublinhado
def titulo(str_titulo):
    print('\n')
    print(str_titulo)
    print(len(str_titulo) * '-')
    
    
# %% Função main - programa principal
def main():
    # %% Apresentação do programa

    titulo('TSL - Teoria Secular Linear')
    
    apresentacao =""" O objetivo deste programa é aplicar a Teoria Secular Linear a
     um exemplo extraído de 'Carpino Teoria lineare delle perturbazioni secolari'
     (o texto está em italiano). O arquivo pode ser obtido no endereço:
     http://www.brera.mi.astro.it/~carpino/didattica/lagrange.pdf."""
     
    print(apresentacao)
    
    # %% Constantes
    # Constante da gravitação universal
    G = (0.01720209895)**2 # para o SS AU^3 d^-2 M_sol^-1
    
    # Massa do Sol
    M = 1.00000598 # Massa 'unitária' para o Sol, segundo Carpino
    
    # %% Dados dos planetas
    
    # Leitura dos dados
    pl = pd.read_fwf('planetas-Carpino.txt')
    
    # Semieixo maior
    a = pl['a'] # unidade AU
    
    # Comprimento do vetor de semieixos.
    len_a = len(a)
    
    # Criando a coluna varpi
    # Se a valor de varpi for superior a 360, subtrair 360 graus.
    pl['varpi'] = pl['capom'] + pl['omega']
    for i in pl.index:
        if pl['varpi'][i] >=360:
            pl.loc[i,'varpi'] = pl.loc[i,'varpi'] - 360
    
        m = pl['m'] # Considerando massa unitária para o Sol
    
    # Dados dos movimentos médios n calculado
    n = (G * (M + m)/a**3)**(1/2) # unidade rad/day
    
    # Dados dos movimentos médios lidos
    #n = pl['n'] # deg/day
    #n = (n * np.pi / 180) # unidade rad/day
    
    # Excentricidade
    e = pl['e']
    
    # Inclinação
    inc = pl['inc'] # unidade deg
    
    # Longitude do periastro
    varpi = pl['varpi'] # unidade deg
    
    # Longitude do nodo ascendente
    capom = pl['capom'] # unidade
    
    
    # %% Verificando dados de entrada
    
    titulo('Dados de entrada')
    print(pl)
    
    titulo('Observação: Unidades')
    print('    [m]: unidade de massa solar')
    print('    [a]: AU')
    print('  [inc]: deg')
    print('[capom]: deg')
    print('[omega]: deg')
    print('    [n]: rad/day')
    print('[varpi]: deg')
    
    titulo('Observação: Data')
    print('Os ângulos foram tomados na data 28 de Julho de 1969,')
    print('no site do JPL da NASA, para ficar de acordo com')
    print('o texto de Carpino, na Tabela III, página 27.')
    
    # %% Verificando função alpha
    # A resultado pode ser verificado no exemplo de Carpino
    # (http://www.brera.mi.astro.it/~carpino/didattica/lagrange.pdf),
    # página 25. A diagonal não está definida e portanto foi 
    # atribuido o valor zero. Como era de se esperar, a matriz é diagonal.
    
    M_alpha = np.zeros((len_a, len_a))
    for i in range(0, len_a):
        for j in range(0,len_a):
            if j != i:
                M_alpha[i, j] = al(a[i], a[j])
            else:
                M_alpha[i, j] = 0
    
    titulo('Função alfa (página 25)')        
    print(M_alpha)  
    
    # %% Verificando função para os Coeficientes de Laplace 
    # b_{3/2}^{1}(\alpha_{i,j})$
    
    M_laplace = np.zeros((len_a, len_a))
    for i in range(0, len_a):
        for j in range(0,len_a):
            if j != i:
                M_laplace[i, j] = coeff_laplace(1, al(a[i],a[j]))
            else:
                M_laplace[i, j] = 0
    titulo('Coeficientes de Laplace para k = 1 (página 25)')
    print(M_laplace)
    
    # %% Verificando função para os Coeficientes de Laplace $b_{3/2}^{2}(\alpha_{i,j})$
    # A resultado pode ser verificado no exemplo de Carpino
    # (http://www.brera.mi.astro.it/~carpino/didattica/lagrange.pdf),
    # página 26. A diagonal não está definida e portanto foi atribuido o
    # valor zero. Como era de se esperar, a matriz é diagonal.
    
    M_laplace = np.zeros((len_a, len_a))
    for i in range(0, len_a):
        for j in range(0,len_a):
            if j != i:
                M_laplace[i, j] = coeff_laplace(2, al(a[i],a[j]))
            else:
                M_laplace[i, j] = 0
                
    titulo('Coeficientes de Laplace para k = 2 (página 26)')
    print(M_laplace)
    
    # %% Determinando a matriz A
    A = np.zeros((len_a, len_a))
    for i in range(0, len_a):
        for j in range(0,len_a):
            if j != i:
                A[i, j] = - parte1(i, G, n, a) * parte2(i, j, m, al, a) * \
                parte3(2, i, j, coeff_laplace, a)
            else:
                parcial = 0
                for l in range(0, len_a):
                    if l != i:
                        parcial = parcial + parte2(i, l, m, al, a) * \
                        parte3(1, i, l, coeff_laplace, a)
                A[i, j] = parte1(i, G, n, a) * parcial
    
    
    # Convertendo a matriz A para unidades de arcsec / yr
    A = A * (180/np.pi) * 365.25 * 3600
    
    # Convertendo a matriz A para unidades de deg / yr
    #A = A * (180/np.pi) * 365.25
    
    # Obtendo a matriz simétrica A
    A_cal = np.zeros((len_a,len_a))
    for i in range(0,len_a):
        for j in range(0,len_a):
            if i == j:
                A_cal[i, i] = A[i,i]
            else:
                A_cal[i, j] = ((a[i] * np.sqrt(m[i] * n[i]))                            / (a[j] * np.sqrt(m[j] * n[j]))) * A[i, j]
                
    # Obtendo os autovetores e autovalres de A
    A_eigen = np.linalg.eig(A_cal)
    
    
    # %% Resultados da matriz A
    # A matriz simétrica, na página 26 do texto de Carpino.
    
    titulo('Resultado da matriz A em arcsec/yr')
    
    print('Matriz A')
    print(A)
    print('\n Matriz A simétrica, conforme Carpino, página 26')
    print(A_cal)
    print('\n Autovalores de A simétrica')
    print(A_eigen[0])
    print('\n Autovetores de A simétrica')
    print(A_eigen[1])
    
    # %% Matriz B
    B = np.zeros((len_a, len_a))
    for i in range(0, len_a):
        for j in range(0,len_a):
            if j != i:
                B[i, j] =  parte1(i, G, n, a) * parte2(i, j, m, al, a) * \
                parte3(2, i, j, coeff_laplace, a)
            else:
                B[i, i] = -A[i, i]
    
    # Convertendo a matriz B para unidades de arcsec / yr
    B = B * (180/np.pi) * 365.25 * 3600
    
    # Convertendo a matriz B para unidades de deg / yr
    #B = B * (180/np.pi) * 365.25
    
    # Obtendo a matriz simétrica B
    B_cal = np.zeros((len_a,len_a))
    for i in range(0,len_a):
        for j in range(0,len_a):
            if i == j:
                B_cal[i, i] = -A[i,i]
            else:
                B_cal[i, j] = ((a[i] * np.sqrt(m[i] * n[i]))                            / (a[j] * np.sqrt(m[j] * n[j]))) * B[i, j]
    
    # Obtendo os autovetores e autovalres de B
    B_eigen = np.linalg.eig(B_cal)
    
    # %% Resultados da matriz B
    # A matriz simétrica, na página 26 do texto de Carpino.
    
    titulo('Resultado da matriz B em arcsec/yr')
    print('Matriz B')
    print(B) 
    print('\n Matriz B simétrica, conforme Carpino, página 26')
    print(B_cal)
    print('\n Autovalores de B simétrica')
    print(B_eigen[0])
    print('\n Autovetores de B simétrica')
    print(B_eigen[1])
    
    # %% Matrizes rescalonadas
    
    # Ainda está muito ruim
    u_estrela = A_eigen[1]
    u = np.zeros((len_a, len_a))
    
    for i in range(0, len_a):
        for j in range(0, len_a):
            u[i, j] = (u_estrela[i,j]/(a[i] * np.sqrt(m[i] * n[i])))
    
    # print(u)
    
    # %% Criando novas colunas com as variáveis regularizadas
    pl['h'] = e * np.sin(np.radians(varpi))
    pl['k'] = e * np.cos(np.radians(varpi))
    pl['p'] = np.sin(np.radians(inc)) * np.sin(np.radians(capom))
    pl['q'] = np.sin(np.radians(inc)) * np.cos(np.radians(capom))
        
    titulo('Variáveis regulares')
    print(pl[['Planeta', 'h', 'k', 'p', 'q']])
    
# %% Execução da função main
if __name__ == '__main__':
    main()

