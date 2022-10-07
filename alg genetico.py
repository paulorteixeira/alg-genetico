import numpy as np
import random as rd
import matplotlib.pyplot as plt

tamCromossomo = 30 ##tamanho do cromossomo
pc = 0.95          ##prop de cruzamento
pm = 0.1           ##prop de mutações
numgeracoes  = 100 ## numero de gerações
tamPopulacao = 50 ## tamanho da populacao
roleta_torneio = False
pontosCruzamento = 2
kpop = 15
elitismo = False

## geração aleatoria da populaçao inicial

p = np.zeros((tamPopulacao,tamCromossomo))##populaçao
for i in range(tamPopulacao):
    for j in range (tamCromossomo):
        a = rd.uniform(0,1)
        if(a>=0.5):
            p[i][j] = 1
        else:
            p[i][j] = 0

## Criação de vars do AG
ind = np.zeros(tamCromossomo)  
individuo = np.zeros(tamPopulacao) ##valores reais
Aptidao = np.zeros(tamPopulacao)
novageracao = np.zeros((tamPopulacao,tamCromossomo))
geracoes = 0

## iniciando o AG
while (geracoes<=numgeracoes):
    novosindividuos = 0
    while (novosindividuos<(tamPopulacao-1)):
        ##Tranformando individuos de bin para real
        for i in range(tamPopulacao):
            ind[:] = p[i,:]
            conv = 0
            for j in range(tamCromossomo):
                conv=conv+ind[j]*(2**(tamCromossomo-(j+1)))## convertendo a base 2 para 10
            individuo[i] = (512/(2**tamCromossomo-1))*conv
        ##calculo da aptidao dos individuos
        TotalAptidao = 0
        for i in range(tamPopulacao):
            Aptidao[i] = abs(individuo[i]*np.sin(np.sqrt(abs(individuo[i]))))+5 ## adaptação da função f(x)
            TotalAptidao = Aptidao[i]+TotalAptidao
        ## Selecao dos pais para cruzamento  
        ## identificar a probabilidade de cada individuo 
        pic = np.zeros(tamPopulacao)   
        pitotal = np.zeros(tamPopulacao)
        pic = (1/TotalAptidao)*Aptidao
        #criando roleta
        if(roleta_torneio == True):
            for i in range(tamPopulacao):
                if(i == 0):
                    pitotal[i] = pic[i]
                else:
                    pitotal[i] = pic[i]+pitotal[i-1]
            ## sorteando os pais de acordo com a prob
            roleta1 = rd.uniform(0,1)
            i = 0
            while (roleta1>pitotal[i]):
                i = i+1
            pai1 = i

            roleta2 = rd.uniform(0,1)
            i = 0
            while (roleta2>pitotal[i]):
                i = i+1
            pai2 = i

            while(pai2 == pai1):
                roleta2 = rd.uniform(0,1)
                i = 0
                while (roleta2>pitotal[i]):
                    i = i+1
                pai2 = i
        ##fim roleta
        else:
            auxii = 0
            maxApt1 = 0
            maxApt2 = 0
            k_sorteado = np.zeros(kpop)
            i = 0
            while(i< kpop): 
                k_sorteado[i] = int(round(1+(tamCromossomo-2)*rd.uniform(0,1)))
                i = i +1
           
            i=0
            while(i< kpop):
                
                opi = int(k_sorteado[i])
                if(maxApt1<Aptidao[opi]):
                    maxApt1 = Aptidao[opi]
                    pai1 = int(k_sorteado[i])
                i = i +1
                
            i = 0
            while(i< kpop):
                opi = int(k_sorteado[i])
                if(maxApt2<Aptidao[opi] and maxApt2!=maxApt1):
                    maxApt2 = Aptidao[opi]
                    pai2 = int(k_sorteado[i])
                i = i +1  
            
        ##operação de cruzamento com 1 ponto
        if(pontosCruzamento == 1):
            if(pc>rd.uniform(0,1)):
                c = round(1+(tamCromossomo-2)*rd.uniform(0,1))
                gene11 = p[pai1][0:c]
                gene12 = p[pai1][c:tamCromossomo]
                gene21 = p[pai2][0:c]
                gene22 = p[pai2][c:tamCromossomo]
                filho1 = np.concatenate((gene11,gene22),axis=None)
                filho2 = np.concatenate((gene21,gene12),axis=None)

                novageracao[novosindividuos,:] = filho1
                
                if(elitismo == True):
                    i=0
                    maxApto = 0
                    nextPai = 0
                    while(i< kpop):    
                        opi = int(k_sorteado[i])
                        if(maxApto<Aptidao[opi]):
                            maxApto = Aptidao[opi]
                            nextPai = int(k_sorteado[i])
                        i = i +1
                    
                    novageracao[novosindividuos,:] = p[nextPai][:]

                novosindividuos = novosindividuos+1    
                novageracao[novosindividuos,:] = filho2                        
                novosindividuos = novosindividuos+1
        else:
            if(pc>rd.uniform(0,1)):
                
                c = round(1+(tamCromossomo-2)*rd.uniform(0,1))## cruzamento entre 2 pontos
                c2 = round(1+(tamCromossomo-2-c)*rd.uniform(0,1)) + c
                
                gene11 = p[pai1][0:c]
                gene12 = p[pai1][c:c2]
                gene13 = p[pai1][c2:tamCromossomo]
                
                gene21 = p[pai2][0:c]
                gene22 = p[pai2][c:c2]
                gene23 = p[pai2][c2:tamCromossomo] 
                
                filho1 = np.concatenate((gene11,gene22,gene13),axis=None)
                filho2 = np.concatenate((gene21,gene12,gene23),axis=None)

                novageracao[novosindividuos,:] = filho1
                
                if(elitismo == True):
                    i=0
                    maxApto = 0
                    nextPai = 0
                    while(i< kpop):    
                        opi = int(k_sorteado[i])
                        if(maxApto<Aptidao[opi]):
                            maxApto = Aptidao[opi]
                            nextPai = int(k_sorteado[i])
                        i = i +1
                    
                    novageracao[novosindividuos,:] = p[nextPai][:]
                    
                novosindividuos = novosindividuos+1
                novageracao[novosindividuos,:] = filho2
                novosindividuos = novosindividuos+1


        
            
            
        ##mutação
        if(pm>rd.uniform(0,1)):
            d = round(1+(tamCromossomo-2)*rd.uniform(0,1))
            if(novageracao[novosindividuos-2][d]==0):
                novageracao[novosindividuos-2][d]=1
            else:
                novageracao[novosindividuos-2][d]=0
            if(novageracao[novosindividuos-1][d]==0):
                novageracao[novosindividuos-1][d]=1
            else:
                novageracao[novosindividuos-1][d]=0

    indice = Aptidao.argmax()
    elem = individuo[indice]

    
    if(geracoes>0):
        if(elem>melhoraptidao):
            melhoraptidao = elem
            indicemelhoraptidao = indice
            melhorx = individuo[indice]
            y = -abs(melhorx*np.sin(np.sqrt(abs(melhorx))))
            melhorgeracao = geracoes
            populacao = novageracao
    else:
        melhoraptidao = elem
        indicemelhoraptidao = indice
        melhorx  = individuo[indice]
        y = -abs(melhorx*np.sin(np.sqrt(abs(melhorx))))
        melhorgeracao = 0
        populacao = p
    
        
    p = novageracao                        
    geracoes = geracoes + 1

    
#gerando o grafico //inf- sup-passos
vetorx = np.linspace(0,450,100)
vetory = -abs(vetorx*np.sin(np.sqrt(abs(vetorx))))
plt.plot(vetorx,vetory)
plt.plot(melhorx,y,'o',color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

print('melhor geracao:',melhorgeracao)
print('melhor x:', melhorx)
print('y:',y)






            
