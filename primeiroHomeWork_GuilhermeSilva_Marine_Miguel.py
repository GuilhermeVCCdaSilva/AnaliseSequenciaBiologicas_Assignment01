#! /usr/bin/python
import sys
from Bio import Entrez


def obterArgumentosDoUtilizador(db,term):
    ''' 

    Obtem os argumentos "Data Base" e o "term" escritos pelo utilizador
    return: 
    Os argumentos "Data Base" db (String) e o "term" term (String) escritos pelo utilizador

    '''

    db = sys.argv[1]
    term = sys.argv[2]
    return db, term


def obterResultadoESearch(db, term):
    ''' 

    Obtem o resultado de um ESearch feito à "Entrez API" com a "Data Base" e o "term" escritos pelo utilizador
    param: 
    db ,term: "Data Base" e o "term" escritos pelo utilizador
    return: 
    o resultado da pesquisa feito ESearch feito à "Entrez API"

    '''

    eSearch = Entrez.esearch(db=db, term=term, usehistory="y")
    resultado = Entrez.read(eSearch)
    return resultado


def obterWebEnv(resultado):
    ''' 
    
    Obtem o WebEnv a partir do resultado da função anterior
    param: 
    resultado: o resultado da pesquisa feito ESearch feito à "Entrez API"
    return: 
    O WebEnv do resultado

    '''

    webEnv = resultado["WebEnv"]
    return webEnv


def obterQueryKey(resultado):
    ''' 
    
    Obtem a QueryKey a partir do resultado da primeira função
    param: 
    resultado: o resultado da pesquisa feito ESearch feito à "Entrez API"
    return: 
    A QueryKey do resultado

    '''

    queryKey = resultado["QueryKey"]
    return queryKey


def obterInformacaoEFetch(queryKey, webEnv):
    ''' 
    
    Obtem a informação do EFetch a partir da "Data Base", da QueryKey e do WebEnv que veem das funções anteriores
    param: 
    queryKey, webEnv: a QueryKey e o WebEnv do resultado feito com o ESearch feito à "Entrez API"
    return: 
    A informacao pelo feito EFetch feiro com a "Data Base", QueryKey e WebEnv no formato fasta 

    '''

    fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey, rettype='fasta')
    informacao = fetchHandle.read()
    return informacao


def escreverFicheiroFasta(informacao):
    ''' 

    Escreve a informação da função anterior num ficheiro Fasta na pasta onde se encontra o programa
    param: 
    A informacao pelo feito EFetch feiro com a "Data Base", QueryKey e WebEnv no formato fasta

    '''
    saveFasta = open(r'sequenciacao.fasta', 'w+')
    saveFasta.write(informacao)
    saveFasta.close()


if __name__ == '__main__':
    db = ""
    term = ""
    db,term = obterArgumentosDoUtilizador(db,term)
    resultado = obterResultadoESearch(db, term)
    queryKey, webEnv = obterQueryKey(resultado), obterWebEnv(resultado)
    informacao = obterInformacaoEFetch(queryKey, webEnv)
    print(informacao)
    #escreverFicheiroFasta(informacao)
    #sys.stderr.write("Foi adicionado com sucesso um ficheiro sequenciacao.Fasta com o output >> ")
    
    #Feito por:
    #Marine Fournier 202000224
    #Guilherme Silva 202000178
    #Miguel 202101030
    
