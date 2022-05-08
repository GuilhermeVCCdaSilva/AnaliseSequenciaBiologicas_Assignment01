#!/usr/bin/bash
escreverSequenciasCoelhos(){

    '''
    Escreve as todas sequências dos Coelhos a partir do primeiroHomeWork_GuilhermeSilva_Marine_Miguel.py 
    com a base de dados (bd) e os (termos) terms como argumentos do programa python.

    params:
    Aceita como parâmetros o primeiro termo, os restantes termos e o nome do ficheiro a ser escrito.
    '''

    bd="Nucleotide"
    primeiroTermo=$1
    terms=$2
    ficheiro=$3
    python3 primeiroHomeWork_GuilhermeSilva_Marine_Miguel.py $bd $primeiroTermo > $ficheiro
    clear
    for i in "${terms[@]}"; do
        python3 primeiroHomeWork_GuilhermeSilva_Marine_Miguel.py $bd $i >> $ficheiro
        clear
    done    
}

escreverSequenciasCoelhosJapCytb(){

   '''
    Executa a função Sequencias escreverSequenciasCoelho() com o primeiroTermo, terms e nome do ficheiro para as sequências
    dos Coelhos do japão com o gene Cytb.
    '''

    primeiroTermo="AB058606.1"
    terms=("AB058605.1" "AB058604.1" "AB058616.2" "AB058615.2" "AB058607.1")
    ficheiro=japCoelhosCytb.fasta
    escreverSequenciasCoelhos $primeiroTermo $terms $ficheiro
}

escreverSequenciasCoelhosJap12S(){
    
    '''
    Executa a função Sequencias escreverSequenciasCoelho() com o primeiroTermo, terms e nome do ficheiro para as sequências
    dos Coelhos do japão com o gene 12S.
    '''
    primeiroTermo="AB058609.1" 
    terms=("AB058608.1" "AB058614.1" "AB058612.1" "AB058613.1" "AB058610.1" "AB058611.1")
    ficheiro=japCoelhos12S.fasta
    escreverSequenciasCoelhos $primeiroTermo $terms $ficheiro
  
}

escreverSequenciasCoelhosWorldWide12S(){

    '''
    Executa a função Sequencias escreverSequenciasCoelho() com o primeiroTermo, terms e nome do ficheiro para as sequências
    dos Coelhos ao redor do mundo com o gene 12S.
    '''

    primeiroTermo="AF176588"
    terms=("U59264" "U31044" "AB053205" "AY012126" "AB053258" "U58921" "U58922")
    ficheiro=coelhosWorldWide12S.fasta
    escreverSequenciasCoelhos $primeiroTermo $terms $ficheiro
}

escreverSequenciasCoelhosWorldWideCytb(){

    '''
    Executa a função Sequencias escreverSequenciasCoelho() com o primeiroTermo, terms e nome do ficheiro para as sequências
    dos Coelhos ao redor do mundo com o gene Cytb.
    '''

    primeiroTermo="AF010156"
    terms=("AF010152" "AF010153" "U58933.1" "AF010161" "AF010154" "AF010155" "AF009732" "AF009733" "U07566" "U58935" "AB053206" "AB053257" "U58930" "U58931")
    ficheiro=coelhosWorldWideCytb.fasta
    escreverSequenciasCoelhos $primeiroTermo $terms $ficheiro

}

escreverFicheiroAlinhamentoCytb(){
    '''
    Faz-se um conjunto de sed (substituições) nos ficheiros coelhosWorldWideCytb.fasta e japCoelhosCytb.fasta
    que são guardos em dois ficheiros temporários. 

    Os ficheiros temporarários juntam-se num só ficheiro temporário final. 
    
    Em seguida é executado um alinhamento (mafft) das sequências do ficheiro temporário final
    e o alinhamento é guaradado num ficheiro denominado alinCytb.fasta.

    No final todo os ficheiros temporários são removidos.
    '''

    sed 's/'mito'.*//' coelhosWorldWideCytb.fasta | sed 's/'cyto'.*//' | sed 's/.*'.1'//' | sed 's/ />/' | sed 's/ /_/' | sed '0,/'timidus'/s//timid_Russia/' | sed '0,/'timidus'/s//timid_UK/' > tempCytbWW.fasta
    sed 's/.*:/?/' japCoelhosCytb.fasta | sed 's/?/>/' | sed 's/ /_/'  > tempCytbJap.fasta
    cat tempCytbWW.fasta >> tempCytbJap.fasta
    mafft tempCytbJap.fasta > alinCytb.fasta
    rm tempCytbWW.fasta
    rm tempCytbJap.fasta
}

escreverFicheiroAlinhamento12S(){
    '''
    Faz-se um conjunto de sed (substituições) nos ficheiros coelhosWorldWide12S.fasta e japCoelhos12S.fasta
    que são guardos em dois ficheiros temporários. 

    Os ficheiros temporarários juntam-se num só ficheiro temporário final. 
    
    Em seguida é executado um alinhamento (mafft) das sequências do ficheiro temporário final
    e o alinhamento é guaradado num ficheiro denominado alin12S.fasta.

    No final todo os ficheiros temporários são removidos.
    '''

    sed 's/'mito'.*//' coelhosWorldWide12S.fasta | sed 's/'12S'.*//' | sed 's/.*'.1'//' | sed 's/ />/' | sed 's/','.*//' | sed 's/'R.680'.*//' | sed 's/ /_/' > temp12SWW.fasta
    sed 's/.*:/?/' japCoelhos12S.fasta | sed 's/?/>/' | sed 's/ /_/' > temp12SJap.fasta
    cat temp12SWW.fasta >> temp12SJap.fasta
    mafft temp12SJap.fasta > alin12S.fasta 
    rm temp12SWW.fasta
    rm temp12SJap.fasta


}

escreverAlinhamentoNexus12S(){

    '''
    Converte o ficheiro denominado alin12S.fasta para nexus: 12S.nex.
    '''
    clear
    seqmagick convert --output-format nexus --alphabet dna alin12S.fasta 12S.nex
    clear
}

escreverAlinhamentoNexusCytb(){ 
    '''
    Converte o ficheiro denominado alin12S.fasta para nexus: Cytb.nex.
    '''
    clear
    seqmagick convert --output-format nexus --alphabet dna alinCytb.fasta Cytb.nex
    clear
}

calcularArvores(){
    '''
    *PAUP permite inferir árvores evolutivas é um programa de filogenia computacional.
    As arvores dos ficheiros nex são calculadas em função aos métodos referidos em cada ficheiro de confuguração.
    E assim os ficheiros bootstrap_MP12S.tre, bootstrap_NJ12S.tre e bootstrap_NJCytb.tre com as informações do melhor 
    bootstrap são criados para as três árvores.
    '''

    paup configMP12SPAUP.txt
    paup configNJ12SPAUP.txt
    paup configNJCytbPAUP.txt
}


escreverArvoreNHKeDesenharToyTree(){
    '''
    O grep permite pesquisar a linha que começa com tree. 
    A partir do sed essa linha é escrita no formato newick.
    Essa escrita é guardada num ficheiro .nhk.
    Isto acontece nas três árvores.
    Por fim é executado o programa toyTree.py que desenha as árvores e guarda em .pdf.
    '''

    grep "^tree" bootstrap_MP12S.tre | sed 's/.*'U]'//' | sed 's/ //' > bootstrap_MP12S.nhk
    grep "^tree" bootstrap_NJ12S.tre| sed 's/.*'U]'//' | sed 's/ //' > bootstrap_NJ12S.nhk
    grep "^tree" bootstrap_NJCytb.tre | sed 's/.*'U]'//' | sed 's/ //' > bootstrap_NJCytb.nhk
    python3 toyTree.py 

}

organizarFicheiros(){
    '''
    Os ficheiros criados são organizados em pastas por meio de mkdir (criação dos ficheiros) e mv (movidos para as pastas).
    Depois os ficheiros e as pastas são apagados (rm) para não existir nenhuma duplicação ou erros em execusões em simultâneo.
    '''

    rm -r originSeq
    mkdir -p originSeq
    rm -r alinSeq
    mkdir -p alinSeq
    rm -r alinNexSeq
    mkdir -p alinNexSeq
    rm -r bootstrapTrees
    mkdir -p bootstrapTrees
    rm -r pdfTrees
    mkdir -p pdfTrees
    mv -iv japCoelhosCytb.fasta originSeq 
    rm japCoelhosCytb.fasta 
    mv -iv japCoelhos12S.fasta originSeq
    rm japCoelhos12S.fasta
    mv -iv coelhosWorldWide12S.fasta originSeq
    rm coelhosWorldWide12S.fasta
    mv -iv coelhosWorldWideCytb.fasta originSeq
    rm coelhosWorldWideCytb.fasta
    mv -iv alin12S.fasta alinSeq
    rm alin12S.fasta
    mv -iv alinCytb.fasta alinSeq
    rm alinCytb.fasta
    mv -iv 12S.nex alinNexSeq
    rm 12S.nex
    mv -iv Cytb.nex alinNexSeq
    rm Cytb.nex
    mv -iv bootstrap_MP12S.nhk bootstrapTrees
    rm bootstrap_MP12S.nhk
    mv -iv bootstrap_MP12S.tre bootstrapTrees
    rm bootstrap_MP12S.tre
    mv -iv bootstrap_NJ12S.nhk bootstrapTrees
    rm bootstrap_NJ12S.nhk
    mv -iv bootstrap_NJCytb.nhk bootstrapTrees
    rm bootstrap_NJCytb.nhk
    mv -iv bootstrap_NJCytb.tre bootstrapTrees
    rm bootstrap_NJCytb.tre
    mv -iv bootstrap_NJ12S.tre bootstrapTrees
    rm bootstrap_NJ12S.tre
    mv -iv tree-plot_bootstrap_NJ12S.pdf pdfTrees
    rm tree-plot_bootstrap_NJ12S.pdf
    mv -iv tree-plot_bootstrap_NJCytb.pdf pdfTrees
    rm  tree-plot_bootstrap_NJCytb.pdf
    mv -iv tree-plot_bootstrap_MP12S.pdf pdfTrees
    rm tree-plot_bootstrap_MP12S.pdf
    clear

}

escreverSequenciasCoelhosJapCytb
escreverSequenciasCoelhosJap12S
escreverSequenciasCoelhosWorldWide12S
escreverSequenciasCoelhosWorldWideCytb
escreverFicheiroAlinhamentoCytb
escreverFicheiroAlinhamento12S
escreverAlinhamentoNexus12S
escreverAlinhamentoNexusCytb
calcularArvores
escreverArvoreNHKeDesenharToyTree
organizarFicheiros
