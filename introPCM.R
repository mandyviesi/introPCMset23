## Introducao aos metodos comparativos filogeneticos ##
## Rladies Ribeirao Preto - 26/09/2023 ##
## Por: Amanda Vieira da Silva ##
## Contato: @MandyViesi (twitter, github)
## Site: amandaviesi.weebly.com ##


# Pacotes ----------------------------------------------------------------

# PCM
library(ape)
library(picante)
library(geiger)
library(phytools)
library(diversitree)

# Mega trees
library(U.PhyloMaker)
# install.packages("U.PhyloMaker")
# devtools::install_github("jinyizju/U.PhyloMaker", force = TRUE)
library(V.PhyloMaker2)

# Manipulacao de dados
library(dplyr)

# Modelos estatisticos
library(nlme)

# Graficos
library(ggplot2)
library(ggtree)


# 1. Apresentacao -------------------------------------------------------

# 2. Inferencia filogenetica x metodos comparativos filogeneticos (PCM) -------

# Inferencia filogenetica:
## https://www.researchgate.net/publication/322701351_Metodologia_da_Inferencia_Filogenetica
## https://www.youtube.com/watch?v=_btf3gqwbXs

# Filogenias como hipoteses

# O que sao os PCMs

# Perguntas que os PCMs conseguem responder
## Evolucao de atributos: como um determinado atributo evoluiu
## Correlacao entre atributos: como um atributo afeta a evolucao de outro atributo
## Diversificacao: taxas de extincao e especiacao
## Biogeografia: quais clados se diversificaram em quais ambientes


# 3. Arvores filogeneticas -------------------------------

## O que sao e como interpretar --------------------------

# Formato Newick (formato parentetico)
parentetico <- "((((A, B), C), D), E);"

# Transformando em arvore


# Classe da arvore


# Estrutura de um objeto do tipo phylo no R
# Edges: comprimento de ramo. Representa distancia genetica em filogramas e tempo em arvores ultrametricas
# Nnode: numero de nos
# Tip.label: numero de terminais


# Plotando a arvore


# Parametros do plot (type)


# Parametros do plot (edge.color)


#Parametros do plot (tip.color)



## Simulando uma arvore ----------------------------------

# Processo estocastico: pure-birth tree
arv_sim <- pbtree(b = 0.1, d = 0, n = 10)

## Manipulando uma arvore -------------------------------

# Plot da arvore


# Incluindo numero dos terminais


# Incluindo numero dos terminais


# Extraindo clados


# Plotando a arvore


# Extraindo terminais


# Plotando a arvore


# Posso plotar a arvore no ggplot?
ggplot(arv_sim) +
  geom_tree() +
  theme_tree()

# Onde achar filogenia ---------------------------

## Mega trees ------------------------------------

# Exemplo 1: anfibios

# Carregando a arvore arvore
anfibios_megatree <- read.tree('https://raw.githubusercontent.com/megatrees/amphibian_20221117/main/amphibian_megatree.tre')

# Carregando os dados
anfibios_generos <- read.csv('https://raw.githubusercontent.com/megatrees/amphibian_20221117/main/amphibian_genus_list.csv', sep=",")

# Usando U.PhyloMaker
arvore_anfibio <- U.PhyloMaker::phylo.maker(tree = anfibios_megatree,
                                            gen.list = anfibios_generos,
                                            sp.list = anfibios_generos,
                                            nodes.type = 1, scenario = 3)
arvore_anfibio

# Pegando so a informacao da arvore (arquivo phylo)


# Plot


# Exemplo 2: angiosperma
# Carregando os dados


# Filogenia


# Arvore ultrametrica


# Plot


## Lendo arvore e dados de um arquivo externo --------

# Arvores podem vir em diferentes formatos
# Arquivo txt


# Arquivo phylo


# Arredondando os comprimentos de ramo para duas casas decimais


# E geram resultados iguais



# 4. PCMS ----------------------------------------------

# Origem dos PCMs
## Ausencia de independencia nas amostras

## Contrastes independentes filogeneticos (PIC) --------

# Regressao linear considerando contraste entre duas especies

# Como funciona
## Pega dois terminais. Calcula a diferenca nos valores dos atributos desses terminais
## Remove esses terminais e deixa o ancestral (media ponderada dos dois terminais). Calcula a diferenca pro proximo terminal
## Padronizacao dos contrastes pela variancia

# PIC assume que a quantidade de variação acumula como uma função linear do tempo separando as especies.

# Pior cenario de Felsenstein

# Criando dois clados em formato parentetico
clade_a <- paste(letters[1:13], collapse = ",")
clade_a <- paste(c("(", clade_a, ")"), collapse = "")

clade_b <- paste(letters[14:26], collapse = ",")
clade_b <- paste(c("(", clade_b, ")"), collapse = "")

clades <- paste(c(clade_a, clade_b), collapse = ",")

clades <-paste(c("(", clades, ")"), collapse = "")
clades <- paste(c(clades, ";"), collapse = "")

# Transformando em arvore
felseinstree <- read.tree(text = clades)

# Plot


# Criando comprimento de ramo


# Enraizando a arvore


# Simulando 2 caracteres (selecionando sementes para o codigo ser reproduzivel)


# Selecionando cores para as especies

# Criando um vetor com as especies

# Data frame

# Plot das especies na arvore


# Plot de dispersao dos atributos


# Adicionando as cores por clado


# Regressao desconsiderando a ancestralidade comum


# Resolvendo politomias

# Ajustando PIC para caracter X


# Ajustando PIC para caracter y


# Dataframe


# Regressao PIC: remover o intercepto


# Plot PIC


## Regressao filogenetica --------------------

# Conhecida como PGLS (phylogenetic linear regression)

# Para rodar um PGLS, precisamos que a arvore filogenetica seja convertida em um tipo especial de objeto no R chamado de estrutura de correlação. Essa estrutura de correlação sera usada para definir a distribuição dos residuos do modelo.

### Brownian Motion ----------------------------

# O que e o BM
## Modelo de evolucao de caracteres continuos em que os atributos mudam com o tempo

# O que gera o BM
## Deriva genetica
## Selecao fraca
## Selecao que muda ao longo do tempo

# Arvore simulada


# Incorporando o atributo

# Fenograma


# O que acontece se a gente aumentar o sig2

# Arvore


# Dados quantitativos

# Criando uma coluna com os nomes das especies

# Reordenando os dados de acordo com a filogenia


# Pegando os dados


# Vetor com os nomes das especies


# Checando nomes


# Para construir o modelo de correlação, vamos usar o corBrownian
# IMPORTANTE: quando estamos criando a estrutura de correlação, precisamos especificar a ordem dos taxons nos nossos dados (argumento form). Caso contrário, ele considera que nossos dados estão na mesma ordem da arvore



# PGLS entre SVL (comprimento focinho-cauda) e PCI_limbs (comprimento das pernas)


# Comparar com PIC
# PIC é um caso especial de PGLS: a estrutura de correlação dos resíduos = correlação esperada entre espécies é diretamente proporcional a sua fração de ancestralidade comum desde a raiz

# Calculando PIC para SVL


# Calculando PIC para perna


# Fazendo a regressao: lembrar que a regressao precisa passar pela origem


# Mostrar coeficientes do pgls



### Modificando lambda --------------------------

# Um caso de relaxamento do modelo é via introducao de um parametro adicional: lambda. O lambda é um multiplicador da dos elementos fora da diagonal da matriz.
# Para um modelo de quadrados minimos ordinarios, desconsiderando a filogenia, lambda = 0, para um PGLS, lambda = 1
# Nos podemos estimar o lambda por um procedimento de maxima verossimilhança.


# Filogenia em estrela


# Plot


# Ajustando PGLS
# ML estimou lambda como 1,01



### Ornstein-Uhlenbeck ----------------------

#Simulando uma arvore com 8 terminais


# Adicionando tempo de divergencia


# Atribuindo os caracteres unicos para cada terminal


# Criando os atributos sob evolucao OU
# alpha: parametro de restricao
# Quando alpha = 0, temos um modelo BM



# Criando cladograma


# Aplicando OU nos Anolis



# PGLS com estrutura de correlacao OU


### ANOVA filogenetica ---------------------

# Atributo continuo (variavel resposta)


# Atributo discreto (variavel preditora)


# Anova filogenetica



### ANCOVA filogenetica ---------------------

# Ancova filogenetica


# Intervalos de confianca


## Evolucao de atributos continuos ------------

# Atributo SVL


### Brownian Motion ---------------------------

# Modelo BM


### Ornstein-Uhlenbeck -------------------------

# Modelo OU


### Early burst -------------------------------

# Modelo de radiacao adaptativa
# Maior variacao na raiz da arvore


### Comparacao de modelos ---------------------

# Calculando AIC


# Peso de Akaike


## Evolucao de atributos discretos ------------

# Atributo


### Equal Rates ------------------------------------

# Ajustando modelo ER (equal rates): taxas de transicao sao iguais


### Symmetric --------------------------------------

# Ajustando modelo SYM (symmetric): taxas de transicao sao simetricas


### All rates different ----------------------------

# Ajustando modelo ARD (all rates different): taxas de transicao sao diferentes


### Selecao de modelos ----------------------------

# Calculando AIC


# Peso de Akaike


### Correlacao evolutiva ----------------------

# Caracteres binarios

# Transformando ecomorfo e SVL em caracteres binarios


# SVL binario


# Ecomorfo binario


# Ajustando correlacao de Pagel


# Plot da matriz de transicao


# Criando um dataframe so com os atributos discretos


# Atribuindo 0 e 1 aos estados de caracteres


# Criando um plot dos atributos



## Sinal filogenetico -------------------------

# Tendencia de que especies relacionadas se parecam mais umas com as outras que o esperado ao acaso
# O jeito mais comum de testar sinal filogenetico e comparando com BM buscando entender se as especies se parecem mais ou menos que o esperado baseado em um modelo BM de evolucao ao longo do tempo


### Lambda de Pagel ---------------------------


### K de Blomberg -----------------------------


### D de Fritz --------------------------------

# Para atributos binarios

# Criando um vetor com os nomes das especies


# Calculando D de Fritz


# 5. Materiais --------------------------------

# Arvores filogeneticas
## Passaros> https://birdtree.org/
## Mamiferos: https://vertlife.org/data/mammals/

# Livros
## Revell e Harmon, 2022 - Phylogenetic Comparative Methods in R
## Nunn, 2011 - The comparative approach in evolutionary anthropology and biology

# Blog
## Blog do phytools: http://blog.phytools.org/
## Blog em portugues: https://cafecomr.com/


# Videos
## Disciplina "Evolucao fenotipica e metodos filogeneticos comparativos" ministrada pelo prof. dr. Diogo Provete (https://www.youtube.com/watch?v=HML46MSWDGw&list=PLy2rjqiD2VP7EtBBqovfBYGrsMes8Hp7P)

# Material no R
## https://pedrohbraga.github.io/PhyloCompMethods-in-R-workshop/PhyloCompMethodsMaterial.html
