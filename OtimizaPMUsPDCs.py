## Bibliotecas
from pyxparser.parser import AnaredeParser #Biblioteca para leitura de arquivos do anarede
import numpy as np
import pandas as pd
from pyscipopt import Model, quicksum #Otimizador 1
from scipy.optimize import minimize, Bounds
import itertools

class AlocaPMUsPDCs:
    def __init__(self, caminho_sistema, seed, verbose=False):

        self.verbose = verbose #Flag para exibição de mensagens de progresso de execução
        self.seed = np.random.default_rng(seed) #Definição da semente para reprodução dos resultados

        ## Acesso à topologia do sistema:
        self.parser = AnaredeParser()
        self.sistema = self.parser.parse_file(caminho_sistema)

        ## Acesso aos dados do sistema:
        self.dbar = self.sistema.get('DBAR', []) #Todas as barras e suas informações
        self.dlin = self.sistema.get('DLIN', []) #Todas as linhas e suas informações

        ## Numeração (identificação) dos barramentos do sistema:
        self.buses_names = [int(bus['number']) for bus in self.dbar] #Nome dos barramentos do sistema
        self.buses_idx = {bus: idx for idx, bus in enumerate(self.buses_names)} #Dicionário dos Índices (ordem numérica) dos barramentos do sistema
        self.nbus = len(self.buses_names) #Total de barramentos no sistema

        ## Prints dos barramentos:
        print(f"\nO número de barramentos do sistema é de {self.nbus}\n")
        print(f"Os barramentos do sistema são dados por:\n{self.buses_names}\n")
        print(f"Os índices dos barramentos do sistema são dados por:\n{self.buses_idx}\n")

        ## Determinação das vizinhanças entre os barramentos do sistema:
        self.viz = [set() for _ in range(self.nbus)] #Lista de conjuntos vazios, para armazenar os barramentos vizinhos (ÍNDICES) de cada barra alvo
        print(self.viz)
        for line in self.dlin: #Loop em função dos ramos do sistema
            bus_de = int(line["from_bus"]) #Nome da barra de
            bus_para = int(line["to_bus"]) #Nome da barra para
            i, j = self.buses_idx[bus_de], self.buses_idx[bus_para] #Índices das barras de e para
            self.viz[i].add(j) #Adição do índice da barra para à vizinhança da barra de
            self.viz[j].add(i) #Adição do índice da barra de à vizinhança da barra para
        print(self.viz)

        ## Dicionário para organização dos vizinhos de cada barra do sistema (NOMES DAS BARRAS):
        self.dict_buses = {self.buses_names[i]: [self.buses_names[j] for j in sorted(self.viz[i])] for i in range(self.nbus)}
        print(f"Vizinhança:\n{self.dict_buses}\n")

        ## Lista de Arrays para determinar o conjunto de barras oberváveis pela instalação de uma PMU em cada barra i (ÍNDICES). Cada array refere-se a instalação da PMU na barra i correspondente
        self.N = [np.array(sorted(v | {i})) for i, v in enumerate(self.viz)] #Ordenação e união dos conjuntos: barra i com suas vizinhas
        print(f"Índice das vizinhanças observáveis pela instalação de PMUs em cada barra:")
        print([buses_obs.tolist() for buses_obs in self.N])
        print(self.N)

        ## Matriz de cobertura A [NBUS X NBUS] -> (A[i,k] = 1 se a barra i é observada por um medidor localizado na barra k)
        self.A = np.zeros((self.nbus, self.nbus)) #Inicialização por zeros
        for i, n in enumerate(self.N): #Loop em função dos N arrays (cada qual correspondente à PMU em uma barra) 
            self.A[i, n] = 1.0 #Incorporação da observabilidade
        print(f"\nMatriz de cobertura:")
        print(self.A)

        ## Vetor de custos das PMUs (Considera-se custo unitário e igual para as PMUS independente do barramento de alocação)
        self.w = np.ones(self.nbus)

    ## -------------------------------------------- OTIMIZAÇÃO PYSCIPOPT ----------------------------------------------------#

    def resolver_OPP(self, tempo_limite=120, var="B", L=None, N=None, dispositivo="PMUs"): #Função de otimização

        ## Verificação para otimização de PDCs -> Se N=None considera-se a estrutura original do sistema para a otimização de PMUs, se não, há a otimização de PDCs considerando apenas as restrições das barras com PMUs
        if N is None:
            N = self.N

        m = Model("PMU_MINLP") #Criação do ambiente de execução

        ## Prints do processo de otimização:
        if not self.verbose:
            m.hideOutput()
        m.setParam("limits/time", tempo_limite) #Tempo máximo em segundos para obter a solução ótima

        ## Modelagem:

        ## Variáveis (Por padrão as variáveis de decisão são tratadas como binárias):
        x = [m.addVar(vtype=var, lb=0.0, ub=1.0, name=f"x{i}") for i in range(self.nbus)] #Criação das variáveis de decisão associadas às PMUs instaladas em cada barramento
        if L is not None: #Se há limites de canais de comunicação:
            # print(f"Uso do limite de canais de comunicação no processo de otimização")
            y = {(k, j): m.addVar(vtype=var, lb=0.0, ub=1.0, name=f"y{k}_{j}") for k in range(self.nbus) for j in self.viz[k]} #Criação das variáveis de decisão associadas aos barramnetos observados pelas PMUs, exemplo: se y[k,j] = 1, a PMU presente em k observa j
        z = m.addVar(vtype="C", lb=0, name="FOB") #Criação de variável contínua para armazenar o valor da FOB

        ## FOB: X^T * W * X = SUM (w_i * x_i^2) -> Minimização do custo de alocação de PMUs:
        m.setObjective(z, "minimize") #Definição do problema de minimização
        m.addCons(quicksum(self.w[i] * x[i] * x[i] for i in range(self.nbus)) <= z) #Adição da EQ da FOB

        ## Restrições de Igualdade (observabilidade): 
        if L is None: #Se não há limite de canais de comunicação:
        ## Restrições: f_i(x) = (1 - x_i) * PROD_{j em viz(i)} (1 - x_j) = 0, para todo i em self.Nbus
            for n in N: #Loop em função das vizinhanças
                prod = 1 #Inicialização do produtório em 1
                for k in n: #Loop em função das barras da vizinhança (ÍNDICES)
                    prod = prod * (1 - x[int(k)]) # Construção da restrição da barra correspondente à vizinhança n
                m.addCons(prod == 0) #Criação da restrição para cada variável x_i
        else: #Se há limites de canais de comunicação:
            ## Restrições: f_i(x) = (1 - x_i) * PROD_{j em Ω_k} (1-yki) = 0, para todo i em self.Nbus 
            for i in range(self.nbus): #Loop em função dos barramentos
                prod = 1 - x[i] #Primeira parcela da restrição
                for k in self.viz[i]: #Loop em função das barras vizinhas
                    prod = prod * (1 - y[(k, i)]) #Incorporação da parcela do produtório
                m.addCons(prod == 0) #Criação da restrição para cada variável x_i

        ## Novas restrições associadas ao problema com limite de canais de comunicação de corrente (Desigualdade):
        if L is not None: #Se há limites de canais de comunicação:
            # print(f"Uso do limite de canais de comunicação no processo de otimização")
            ## Restrição para impedir que uma barra sem PMU observe suas vizinhas: y_i_j <= x_i, para todo i em self.Nbus e para todo j em Ω_i
            for (k, j), y_var in y.items(): #Loop em função das variáveis y (k indica a barra que pode ter PMU, j a barra vizinhança a barra k correspondente e y_var a variável y para identificar se há observabilidade em j pela presença de PMU em k)
                m.addCons(y_var <= x[k]) #Criação da restrição para cada variável y_i_j
            ## Restrição referente à capacidade de canais de comunicação, de modo que a PMU instalada na barra k pode alocar no máximo L canais dentre os ykj ramos disponíveis: SUM_{j em Ω_k} (y_k_j) = L - x_k, para todo k em self.Nbus
            for k in range(self.nbus): #Loop em função dos barramentos do sistema
                m.addCons(quicksum(y[(k, j)] for j in self.viz[k]) <= L * x[k]) #Criação da restrição para cada variável x_k

        ## Otimização:
        m.optimize()

        ## Condição de não solução:
        if m.getNSols() == 0: 
            print("\nMINLP: nenhuma solução encontrada.")
            return None, None

        ## Para o caso de se encontrar solução ótima:

        ## Recuperação das soluções encontradas pelo SCIP
        solucoes = m.getSols()
        ## Print de todas as soluções encontradas:
        print(f"\nTotal de soluções armazenadas: {len(solucoes)}")
        for i, sol in enumerate(solucoes): #Loop em função das soluções
            x_valores = np.array([m.getSolVal(sol, xi) for xi in x])
            print(f"\nSolução {i+1}:")
            print(f"FOB: {m.getSolVal(sol, z):.10f}")
            print("Valores de x:", np.round(x_valores, 6))
            tol = 1e-5
            barras = [self.buses_names[j] for j in range(self.nbus) if x_valores[j] >= 1 - tol]
            print("Barras com PMUs:", barras)

        ## Apresentação da solução ótima retornada pelo método:
        x_raw = np.array([m.getVal(xi) for xi in x])
        print(f"Desvio máximo da binariedade: {np.max(np.abs(x_raw - np.round(x_raw))):.2e}")
        x_bin = np.round(x_raw).astype(int)
        Fob = m.getVal(z)
        t_total =  m.getTotalTime()
        print(f"Tempo total do solver: {t_total:.20f} s") #Tempo total desde a criação do Model"
        print(f"Nós explorados pelo solver: {m.getNNodes()}") #Nós explorados no branch-and-bound
        barras = self._apresentar(x_bin, Fob, f"MINLP canais L={L} (SCIP, status: {m.getStatus()})", dispositivo)
        return x_bin, Fob, barras

    ## -------------------------------------- SAÍDA (PRINT) SOLUÇÃO ÓTIMA ------------------------------------------- ##

    def _apresentar(self, solucao, fob, titulo, dispositivo):
        if solucao is not None:
            barras = [self.buses_names[i] for i in range(self.nbus) if solucao[i] == 1]
            print(f"\n--- OTIMIZAÇÃO CONCLUÍDA: {titulo} ---")
            print(f"Total MÍNIMO de {dispositivo} alocadas: {int(solucao.sum())} (custo {fob:.2f})")
            print(f"Instalações realizadas nas barras: {barras}")
            if dispositivo == "PMUs":
                print(f"Sistema totalmente observável: {self.observavel(solucao)}")
            return barras
        else:
            print("\nFalha: nenhuma solução viável encontrada.")

    def observavel(self, x_bin): #Verificação de Observabilidade:
        ## Verificação se todas as barras são observáveis para um vetor binário
        return bool(np.all(self.A @ x_bin >= 1))

    ## ------------------------ Alocação ótima de PMU-PDC ----------------------------------- ##

    def Otimiza_PDCs(self, buses_pmu, var="B"):

        print(f"\n---------------  OTIMIZAÇÃO DE PDCS -------------------")

        ## Obtenção dos índices das barras nas quais há a instalação de PMUs:
        buses_pmu_idx = [] #Lista para armazenar os índices das barras com PMU
        for bus in buses_pmu: #Loop em função das barras com PMU instaladas
            bus_idx = self.buses_names.index(bus) #Obtenção do índice do barramento com PMU
            buses_pmu_idx.append(bus_idx) #Adição do índice do barramento com PMU na lista de índices

        ## String para identificar a otimização do dispositivo PDC:
        dispositivo = "PDCs"

        ## Obtenção da vizinhança (restrições) para as barras com PMUs:
        N_pdc = [] #Lista para armazenar as novas restrições -> apenas dos barramentos com PMU
        for i, n in enumerate(self.N): #Loop em função das vizinhanças de cada barramento
            if i in buses_pmu_idx: #Se a vizinhança corresponde a de uma barra com PMU (se a barra i tem PMU)
                N_pdc.append(n) #Adição da vizinhança da barra
        # print(N_pdc)

        ## Solução automática:
        self.resolver_OPP(dispositivo=dispositivo, N=N_pdc, var=var)

