import OtimizaPMUsPDCs as OPP

## Parâmetros Globais:
L=2 #Limite de canais das PMUs
seed=0 #Semente de simulação
var = "C" #Necessário para tratar as variáveis como contínuas
## Sistemas para teste:
## Sistema de 14 Barras usado no artigo:
caminho_sistema = r"C:\Cepel\Anarede\V120001\Exemplos\IEEE14.pwf" #14BUS
## Outros sistemas testados:
# caminho_sistema = r"C:\Cepel\Anarede\V120001\Exemplos\IEEE30.PWF" #30BUS
# caminho_sistema = r"C:\Cepel\Anarede\V120001\Exemplos\IEEE57.pwf" #57BUS
# caminho_sistema = r"C:\Cepel\Anarede\V120001\Exemplos\IEEE118.PWF" #118
# caminho_sistema = r"C:\Cepel\Anarede\V120001\Exemplos\NewEngland.PWF" #Sistema NEW ENGLAND (39BUS)

## Criação do objeto:
Teste = OPP.AlocaPMUsPDCs(caminho_sistema,seed=seed)

## Otimização de PMUs:
# X_bin, fob, barras = Teste.resolver_OPP(L=L, var=var) #Exemplo com limite de canais e variável contínua
X_bin, fob, barras = Teste.resolver_OPP(L=L) #Exemplo com limite de canais e variável binária
# X_bin, fob, barras = Teste.resolver_OPP(var=var) #Exemplo sem limite de canais e variável contínua
# X_bin, fob, barras = Teste.resolver_OPP() #Exemplo sem limite de canais e variável binária

## Otimização de PDCs:
## Se nenhum dos exemplos de alocação for utilizado é necessário rodar a alocação de PMUs antes, para obter as barras de instalação
# barras = [1, 5, 6, 8] #Exemplo de uma alocação específica de PMUs para o 14 barras
# barras = [2,8,10,13] #Exemplo de uma alocação específica de PMUs para o 14 barras
print(f"As barras com instalação de PMUs são dadas por:{barras}")
# Teste.Otimiza_PDCs(barras, var=var) #Otimização de PDCs com variável contínua 
# Teste.Otimiza_PDCs(barras) #Otimização de PDCs como variável binária