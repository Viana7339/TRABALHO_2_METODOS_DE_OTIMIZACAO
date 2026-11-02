Para o funcionamento do programa é necessário determinar o diretório no qual o arquivo correspondente ao PWF do sistema elétrico
se encontra em seu computador, alterando-o na função mainPMUsPDCs.py em caminho_sistema = r"C:\Cepel\Anarede\V120001\Exemplos\IEEE14.pwf" #14BUS

Na função OtimizaPMUsPDCs.py, deixando o parâmetro verbose=True é possível obter vários prints de execução do solver de otimização