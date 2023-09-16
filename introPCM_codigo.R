## Introducao aos metodos comparativos filogeneticos ##
## Rladies Ribeirao Preto - 26/09/2023 ##
## Por: Amanda Vieira da Silva ##
## Contato: @MandyViesi (twitter, github)
## E-mail: mandyvieira15@gmail.com ##
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
arv1 <- read.tree(text = parentetico)

# Classe da arvore
class(arv1)

# Estrutura de um objeto do tipo phylo no R
# Edges: comprimento de ramo. Representa distancia genetica em filogramas e tempo em arvores ultrametricas
# Nnode: numero de nos
# Tip.label: numero de terminais
str(arv1)

# Plotando a arvore
plot(arv1)

# Parametros do plot (type)
plot(arv1, type = "fan")

# Parametros do plot (edge.color)
plot(arv1, edge.color = "red")

#Parametros do plot (tip.color)
plot(arv1, tip.color = "red")


## Simulando uma arvore ----------------------------------

# Processo estocastico: pure-birth tree
arv_sim <- pbtree(b = 0.1, d = 0, n = 10)

## Manipulando uma arvore -------------------------------

# Plot da arvore
plot(arv_sim)

# Incluindo numero dos terminais
tiplabels()

# Incluindo numero dos terminais
nodelabels()

# Extraindo clados
arv2 <- extract.clade(phy = arv_sim, node = 15)

# Plotando a arvore
plot(arv2)

# Extraindo terminais
arv3 <- drop.tip(phy = arv_sim, tip = c(1,2))

# Plotando a arvore
plot(arv3)

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
arvore_anfibio <- U.PhyloMaker::phylo.maker(tree = anfibios_megatree, gen.list = anfibios_generos, sp.list = anfibios_generos, nodes.type = 1, scenario = 3)
arvore_anfibio

# Pegando so a informacao da arvore (arquivo phylo)
arvore_anfibio <- arvore_anfibio[[1]]

# Plot
plot(arvore_anfibio)

# Exemplo 2: angiosperma
# Carregando os dados
plantas_nef <- readxl::read_xlsx("arvore_nef_exemplo.xlsx")

# Filogenia
arvore_plantas_nef <- V.PhyloMaker2::phylo.maker(sp.list = plantas_nef)

# Arvore ultrametrica
ultrametrica <-  compute.brlen(arvore_plantas_nef$scenario.3, method = "Grafen")

# Plot
plot(ultrametrica)

## Lendo arvore e dados de um arquivo externo --------

# Arvores podem vir em diferentes formatos
# Arquivo txt
anolis_arv <- read.tree("anolis_tree.txt")

# Arquivo phylo
anolis_arv2 <- read.tree("anolis.phy")

# Arredondando os comprimentos de ramo para duas casas decimais
anolis_arv2$edge.length <- round(x = anolis_arv2$edge.length, digits = 3)

# E geram resultados iguais
all.equal(anolis_arv, anolis_arv2)


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
plot(felseinstree)

# Criando comprimento de ramo
arvore <- compute.brlen(felseinstree)

# Enraizando a arvore
arvore <- root(arvore, node = 27)

# Simulando 2 caracteres (selecionando sementes para o codigo ser reproduzivel)
set.seed(seed = 87)
carac_x <- fastBM(arvore)
set.seed(seed = 23)
carac_y <- fastBM(arvore)

# Selecionando cores para as especies
cores <- c(rep(x = "#00274c", times = 13), rep(x = "#ffcb05", times = 13))

# Criando um vetor com as especies
especies <- letters

# Data frame
caracteres <- data.frame(carac_x, carac_y, cores, especies)

# Plot das especies na arvore
plotTree(arvore, xlim = c(0, 1.3), ylim = c(1,26), ftype = "off")
points(rep(x = 1,times = 13), y = 1:13, pch = 19, col = "#00274c")
points(rep(x = 1,times = 13), y = 14:26, pch = 19, col = "#ffcb05")

# Plot de dispersao dos atributos
ggplot(caracteres, aes(x = carac_x, y = carac_y)) +
  geom_point() +
  labs(x = "Atributo X", y = "Atributo Y") +
  theme_classic()

# Adicionando as cores por clado
ggplot(caracteres, aes(x = carac_x, y = carac_y)) +
  geom_point(col = cores) +
  labs(x = "Atributo X", y = "Atributo Y") +
  theme_classic()

# Regressao desconsiderando a ancestralidade comum
mod1 <- lm(carac_y ~ carac_x, data = caracteres)
summary(mod1)

# Resolvendo politomias
arvore <- multi2di(arvore)

# Ajustando PIC para caracter X
cx <- setNames(caracteres[,1], rownames(caracteres$especies))
pic_cx <- pic(x = cx, phy = arvore)

# Ajustando PIC para caracter y
cy <- setNames(caracteres[,2], rownames(caracteres$especies))
pic_cy <- pic(x = cy, phy = arvore)

# Dataframe
pics <- data.frame(pic_cx, pic_cy)

# Regressao PIC: remover o intercepto
pic_lm <- lm(pic_cy ~ pic_cx - 1, data = pics)
summary(pic_lm)

# Plot PIC
ggplot(data = pics, aes(x = pic_cx, y = pic_cy)) +
  geom_point() +
  geom_vline(xintercept = 0, col = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0, col = "gray", linetype = "dashed") +
  geom_line(aes(y = predict(pic_lm, pics[2], type = "response")), colour = "red") +
  theme_classic()


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
arv_sim <- pbtree(b = 0.1, d = 0, n = 10)

# Incorporando o atributo
atributo1 <- fastBM(arv_sim, a = 0 , sig2 = 0.1, internal=TRUE)

# Fenograma
phenogram(arv_sim, atributo1, spread.labels=TRUE)

# O que acontece se a gente aumentar o sig2

# Arvore
anolis_arvore <- read.tree("anolis.phy")
anolis_arvore

# Dados quantitativos
anolis_dados <- read.csv("anolisDataAppended.csv", row.names = 1)

# Reordenando anolis_dados de acordo com a filogenia
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
arv1 <- read.tree(text = parentetico)

# Classe da arvore
class(arv1)

# Estrutura de um objeto do tipo phylo no R
# Edges: comprimento de ramo. Representa distancia genetica em filogramas e tempo em arvores ultrametricas
# Nnode: numero de nos
# Tip.label: numero de terminais
str(arv1)

# Plotando a arvore
plot(arv1)

# Parametros do plot (type)
plot(arv1, type = "fan")

# Parametros do plot (edge.color)
plot(arv1, edge.color = "red")

#Parametros do plot (tip.color)
plot(arv1, tip.color = "red")


## Simulando uma arvore ----------------------------------

# Processo estocastico: pure-birth tree
arv_sim <- pbtree(b = 0.1, d = 0, n = 10)

## Manipulando uma arvore -------------------------------

# Plot da arvore
plot(arv_sim)

# Incluindo numero dos terminais
tiplabels()

# Incluindo numero dos terminais
nodelabels()

# Extraindo clados
arv2 <- extract.clade(phy = arv_sim, node = 15)

# Plotando a arvore
plot(arv2)

# Extraindo terminais
arv3 <- drop.tip(phy = arv_sim, tip = c(1,2))

# Plotando a arvore
plot(arv3)

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
arvore_anfibio <- U.PhyloMaker::phylo.maker(tree = anfibios_megatree, gen.list = anfibios_generos, sp.list = anfibios_generos, nodes.type = 1, scenario = 3)
arvore_anfibio

# Pegando so a informacao da arvore (arquivo phylo)
arvore_anfibio <- arvore_anfibio[[1]]

# Plot
plot(arvore_anfibio)

# Exemplo 2: angiosperma
# Carregando os dados
plantas_nef <- readxl::read_xlsx("arvore_nef_exemplo.xlsx")

# Filogenia
arvore_plantas_nef <- V.PhyloMaker2::phylo.maker(sp.list = plantas_nef)

# Arvore ultrametrica
ultrametrica <-  compute.brlen(arvore_plantas_nef$scenario.3, method = "Grafen")

# Plot
plot(ultrametrica)

## Lendo arvore e dados de um arquivo externo --------

# Arvores podem vir em diferentes formatos
# Arquivo txt
anolis_arv <- read.tree("anolis_tree.txt")

# Arquivo phylo
anolis_arv2 <- read.tree("anolis.phy")

# Arredondando os comprimentos de ramo para duas casas decimais
anolis_arv2$edge.length <- round(x = anolis_arv2$edge.length, digits = 3)

# E geram resultados iguais
all.equal(anolis_arv, anolis_arv2)


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
plot(felseinstree)

# Criando comprimento de ramo
arvore <- compute.brlen(felseinstree)

# Enraizando a arvore
arvore <- root(arvore, node = 27)

# Simulando 2 caracteres (selecionando sementes para o codigo ser reproduzivel)
set.seed(seed = 87)
carac_x <- fastBM(arvore)
set.seed(seed = 23)
carac_y <- fastBM(arvore)

# Selecionando cores para as especies
cores <- c(rep(x = "#00274c", times = 13), rep(x = "#ffcb05", times = 13))

# Criando um vetor com as especies
especies <- letters

# Data frame
caracteres <- data.frame(carac_x, carac_y, cores, especies)

# Plot das especies na arvore
plotTree(arvore, xlim = c(0, 1.3), ylim = c(1,26), ftype = "off")
points(rep(x = 1,times = 13), y = 1:13, pch = 19, col = "#00274c")
points(rep(x = 1,times = 13), y = 14:26, pch = 19, col = "#ffcb05")

# Plot de dispersao dos atributos
ggplot(caracteres, aes(x = carac_x, y = carac_y)) +
  geom_point() +
  labs(x = "Atributo X", y = "Atributo Y") +
  theme_classic()

# Adicionando as cores por clado
ggplot(caracteres, aes(x = carac_x, y = carac_y)) +
  geom_point(col = cores) +
  labs(x = "Atributo X", y = "Atributo Y") +
  theme_classic()

# Regressao desconsiderando a ancestralidade comum
mod1 <- lm(carac_y ~ carac_x, data = caracteres)
summary(mod1)

# Resolvendo politomias
arvore <- multi2di(arvore)

# Ajustando PIC para caracter X
cx <- setNames(caracteres[,1], rownames(caracteres$especies))
pic_cx <- pic(x = cx, phy = arvore)

# Ajustando PIC para caracter y
cy <- setNames(caracteres[,2], rownames(caracteres$especies))
pic_cy <- pic(x = cy, phy = arvore)

# Dataframe
pics <- data.frame(pic_cx, pic_cy)

# Regressao PIC: remover o intercepto
pic_lm <- lm(pic_cy ~ pic_cx - 1, data = pics)
summary(pic_lm)

# Plot PIC
ggplot(data = pics, aes(x = pic_cx, y = pic_cy)) +
  geom_point() +
  geom_vline(xintercept = 0, col = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0, col = "gray", linetype = "dashed") +
  geom_line(aes(y = predict(pic_lm, pics[2], type = "response")), colour = "red") +
  theme_classic()


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
arv_sim <- pbtree(b = 0.1, d = 0, n = 10)

# Incorporando o atributo
atributo1 <- fastBM(arv_sim, a = 0 , sig2 = 0.1, internal=TRUE)

# Fenograma
phenogram(arv_sim, atributo1, spread.labels=TRUE)

# O que acontece se a gente aumentar o sig2

# Arvore
anolis_arvore <- read.tree("anolis.phy")
anolis_arvore

# Dados quantitativos
anolis_temp <- read.csv("anolisDataAppended.csv", row.names = 1)

# Criando uma coluna com os nomes das especies
anolis_temp$especies <- rownames(anolis_temp)

# Reordenando os dados de acordo com a filogenia
anolis_ordem <- comparative.data(phy = anolis_arvore, data = anolis_temp, names.col = especies)

# Pegando os dados
anolis_dados <- anolis_ordem$data

# Vetor com os nomes das especies
spp <- rownames(anolis_dados)

# Checando nomes
name.check(anolis_arvore, anolis_dados)

# Para construir o modelo de correlação, vamos usar o corBrownian
# IMPORTANTE: quando estamos criando a estrutura de correlação, precisamos especificar a ordem dos taxons nos nossos dados (argumento form). Caso contrário, ele considera que nossos dados estão na mesma ordem da arvore
corBM <- corBrownian(phy = anolis_arvore, form = ~spp)
names(anolis_dados)

# PGLS entre SVL (comprimento focinho-cauda) e PCI_limbs (comprimento das pernas)
pgls_anolis <- gls(SVL ~ PCI_limbs, data = anolis_dados,
                    correlation = corBM)
pgls_anolis
summary(pgls_anolis)

# Comparar com PIC
# PIC é um caso especial de PGLS: a estrutura de correlação dos resíduos = correlação esperada entre espécies é diretamente proporcional a sua fração de ancestralidade comum desde a raiz
names(anolis_dados)

# Calculando PIC para SVL
svl <- setNames(anolis_dados[,"SVL"], rownames(anolis_dados))
pic_svl <- pic(svl, phy = anolis_arvore)

# Calculando PIC para perna
PCI_limbs <- setNames(anolis_dados[,"PCI_limbs"], rownames(anolis_dados))
pic_limbs <- pic(PCI_limbs, phy = anolis_arvore)

# Fazendo a regressao: lembrar que a regressao precisa passar pela origem
lm(pic_svl ~ pic_limbs -1)

# Mostrar coeficientes do pgls
pgls_anolis$coefficients


### Modificando lambda --------------------------

# Um caso de relaxamento do modelo é via introducao de um parametro adicional: lambda. O lambda é um multiplicador da dos elementos fora da diagonal da matriz.
# Para um modelo de quadrados minimos ordinarios, desconsiderando a filogenia, lambda = 0, para um PGLS, lambda = 1
# Nos podemos estimar o lambda por um procedimento de maxima verossimilhança.
corLambda1 <- corPagel(value = 1, phy = anolis_arvore, form = ~spp)
corLambda1

# Filogenia em estrela
filo_estrela <- rescale(anolis_arvore, model = "lambda", 0)

# Plot
plot(filo_estrela)

# Ajustando PGLS
# ML estimou lambda como 1,01
pgls_lambda <- pgls_anolis <- gls(SVL ~ PCI_limbs, data = anolis_dados,
                                  correlation = corLambda1)
summary(pgls_lambda)


### Ornstein-Uhlenbeck ----------------------

#Simulando uma arvore com 8 terminais
arvore_sim2 <-pbtree(n = 8, scale = 100)

# Adicionando tempo de divergencia
arvore_sim2 <- make.era.map(arvore_sim2, limits = 0:200/2)

# Atribuindo os caracteres unicos para cada terminal
arvore_sim2 <- map.to.singleton(arvore_sim2)

# Criando os atributos sob evolucao OU
# alpha: parametro de restricao
# Quando alpha = 0, temos um modelo BM
atributo2 <- fastBM(arvore_sim2, sig2 = 0.1,
                    alpha = 1, theta = 3,internal=TRUE)

# Criando cladograma
phenogram(arvore_sim2, atributo2, spread.labels=TRUE)

# Aplicando OU nos Anolis
corOU <- corMartins(phy = anolis_arvore, form = ~spp, value = 1)

# PGLS com estrutura de correlacao OU
pgls_OU <- gls(SVL ~ PCI_limbs, data = anolis_dados,
                                  correlation = corOU)
summary(pgls_OU)

### ANOVA filogenetica ---------------------

# Atributo continuo (variavel resposta)
SVL <- setNames(anolis_dados[,"SVL"], rownames(anolis_dados))

# Atributo discreto (variavel preditora)
ecomorfo <- setNames(anolis_dados[,"ecomorph"], rownames(anolis_dados))

# Anova filogenetica
anolis_anova <- phylANOVA(tree = anolis_arvore, x = ecomorfo, y = SVL)
anolis_anova


### ANCOVA filogenetica ---------------------

# Ancova filogenetica
anolis_ancova <- gls(SVL ~ PCI_limbs + ecomorph, data = anolis_dados, correlation = corBM)
anolis_ancova

# Intervalos de confianca
intervals(anolis_ancova)

## Evolucao de atributos continuos ------------

# Atributo
SVL

### Brownian Motion ---------------------------

# Modelo BM
mod_BM <- fitContinuous(phy = anolis_arvore, dat = SVL)
mod_BM

### Ornstein-Uhlenbeck -------------------------

# Modelo OU
mod_OU <- fitContinuous(phy = anolis_arvore, dat = SVL, model ="OU")
mod_OU

### Early burst -------------------------------

# Modelo de radiacao adaptativa
# Maior variacao na raiz da arvore
mod_EB <- fitContinuous(phy = anolis_arvore, dat = SVL, model ="EB")
mod_EB

### Comparacao de modelos ---------------------

# Calculando AIC
aic_continuo <- setNames(c(AIC(mod_BM),
                         AIC(mod_OU),
                         AIC(mod_EB)),
                       c("BM", "OU", "EB"))

aic_continuo

# Peso de Akaike
aic.w(aic_continuo)

## Evolucao de atributos discretos ------------

# Atributo
ecomorfo <- setNames(anolis_dados[,"ecomorph"], rownames(anolis_dados))

### Equal Rates ------------------------------------

# Ajustando modelo ER (equal rates): taxas de transicao sao iguais
mod_ER <- fitDiscrete(phy = anolis_arvore, dat = ecomorfo, model = "ER")

### Symmetric --------------------------------------

# Ajustando modelo SYM (symmetric): taxas de transicao sao simetricas
#mod_SYM <- fitDiscrete(phy = anolis_arvore, dat = ecomorfo, model = "SYM")

### All rates different ----------------------------

# Ajustando modelo ARD (all rates different): taxas de transicao sao diferentes
#mod_ARD <- fitDiscrete(phy = anolis_arvore, dat = ecomorfo, model = "ARD")

### Selecao de modelos ----------------------------

# Calculando AIC
aic_discreto<- setNames(c(AIC(mod_ER),
                           AIC(mod_SYM),
                           AIC(mod_ARD)),
                         c("ER", "SYM", "ARD"))

aic_discreto

# Peso de Akaike
aic.w(aic_discreto)

### Correlacao evolutiva ----------------------

# Caracteres binarios

# Transformando ecomorfo e SVL em caracteres binarios
anolis_dados <- anolis_dados |>
  mutate(ecomorphbi = case_when(
    ecomorph == "TG" ~ "TG",
    ecomorph != "TG" ~ "nao-TG"),
    SVLbi = case_when(
      SVL < 4 ~ "pequeno",
      SVL >= 4 ~ "grande"))

# SVL binario
SVLbi <- setNames(anolis_dados[,"SVLbi"], rownames(anolis_dados))

# Ecomorfo binario
ecomorphbi <- setNames(anolis_dados[,"ecomorphbi"], rownames(anolis_dados))

# Ajustando correlacao de Pagel
pagel <- fitPagel(tree = anolis_arvore, x = SVLbi, y = ecomorphbi)

# Plot da matriz de transicao
plot(pagel)

# Criando um dataframe so com os atributos discretos
names(anolis_dados)

atributos <- anolis_dados |>
  dplyr::select(ecomorphbi, SVLbi)

# Atribuindo 0 e 1 aos estados de caracteres
atributos <- atributos |>
  mutate(ecomorphbi = case_when(
    ecomorphbi == "nao-TG" ~ 1,
    ecomorphbi == "TG" ~ 0),
    SVLbi = case_when(SVLbi == "grande" ~ 1,
    SVLbi == "pequeno" ~ 0))

# Criando um plot dos atributos
trait.plot(anolis_arvore, atributos,
           lab = c("Ecomorfo", "SVL"),
           cex.lab = 0.6, cex.legend = 0.7, w = 0.12,
           str = list(ecomorphbi = c("TG", "not-TG"),
                      SVLbi = c("pequeno", "grande")),
           cols=list(ecomorphbi=c("#588bae", "#00274c"),
                     SVLbi=c("#ffcb05", "#E29000")))


## Sinal filogenetico -------------------------

# Tendencia de que especies relacionadas se parecam mais umas com as outras que o esperado ao acaso
# O jeito mais comum de testar sinal filogenetico e comparando com BM buscando entender se as especies se parecem mais ou menos que o esperado baseado em um modelo BM de evolucao ao longo do tempo


### Lambda de Pagel ---------------------------
phylosig(anolis_arvore, SVL, test = TRUE, method = "lambda")

### K de Blomberg -----------------------------
phylosig(anolis_arvore, SVL, test = TRUE, method = "K")

### D de Fritz --------------------------------
library(caper)

# Para atributos binarios

# Criando um vetor com os nomes das especies
species <- rownames(atributos)

#
comparative.data(phy = anolis_arvore, data = atributos, names.col = species)

# Calculando D de Fritz
phylo.d(data = atributos, phy = anolis_arvore, binvar = SVLbi, names.col = species)

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

