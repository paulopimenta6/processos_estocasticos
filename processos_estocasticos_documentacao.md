# Processos Estocásticos: guia intuitivo e rigoroso — com exemplos em R

**Uma documentação completa para quem está começando: da intuição à implementação.**

*Versão 1.0 · Novembro de 2026 · Escrita a partir de fontes públicas citadas ao longo do texto.*

---

## Sobre este documento

Este guia foi construído para levar você, leitor iniciante, do primeiro contato com a ideia de "aleatoriedade que evolui no tempo" até a implementação dos modelos clássicos em linguagem R. Ele foi escrito com duas exigências simultâneas: **rigor conceitual** (as definições estão corretas e as hipóteses são explicitadas) e **linguagem acessível** (cada conceito vem acompanhado de intuição, exemplo e, sempre que útil, código).

Um aviso prático e honesto: o ambiente onde este texto foi redigido não possui uma instalação de R disponível para execução. Todos os códigos foram escritos para rodar em R (≥ 4.2) com o mínimo de dependências externas — a maior parte funciona apenas com o R base —, mas recomendamos que você execute tudo na sua própria máquina. No fim, há um script pronto (`demo_processos_estocasticos.R`) que reproduz os principais exemplos, além de arquivos CSV com os dados reais usados no estudo de caso climático.

### Para quem é este texto

- Estudantes de graduação em estatística, engenharia, ciência de dados, economia, biologia ou áreas afins que estão vendo processos estocásticos pela primeira vez;
- Profissionais que já usam R no dia a dia e querem entender os modelos probabilísticos por trás de simulações, séries temporais e análises espaciais;
- Qualquer pessoa curiosa sobre como se modela o acaso — de preços de ações a epidemias, de filas de banco a chuvas em São Paulo.

### Pré-requisitos mínimos

Você vai aproveitar mais o texto se tiver familiaridade básica com: probabilidade (eventos, probabilidade condicional, independência), variáveis aleatórias (esperança, variância) e programação elementar em R (vetores, funções, `for`, plots). Se algo faltar, as seções 2 e 3 relembram o essencial sem formalismo excessivo.

### Como usar

- Leia as seções 1–4 como um curso: cada classe de processo segue o roteiro *intuição → definição → propriedades → simulação em R*.
- A seção 5 ensina a parte prática de simulação estocástica em R (geração de números aleatórios, reprodutibilidade, pacotes).
- A seção 6 traz aplicações guiadas completas: finanças, epidemias, filas, clima e geoespacial.
- A seção 7 discute armadilhas, erros comuns e controvérsias — leitura obrigatória para quem quer usar esses métodos com responsabilidade.
- As seções finais mapeiam materiais em português para continuar estudando, referências e limitações deste documento.

### Sumário

1. [A ideia central: o que é um processo estocástico](#1-a-ideia-central-o-que-um-processo-estocstico)
2. [Ferramentas essenciais de probabilidade](#2-ferramentas-essenciais-de-probabilidade)
3. [Informação no tempo: filtrações, tempos de parada e martingalas](#3-informao-no-tempo-filtraes-tempos-de-parada-e-martingalas)
4. [As grandes famílias de processos](#4-as-grandes-famlias-de-processos)
5. [Simulação estocástica em R: do básico ao reproduzível](#5-simulao-estocstica-em-r-do-bsico-ao-reproduzvel)
6. [Aplicações guiadas](#6-aplicaes-guiadas)
7. [Armadilhas, erros comuns e controvérsias](#7-armadilhas-erros-comuns-e-controvrsias)
8. [Trilha de estudos e recursos em português](#8-trilha-de-estudos-e-recursos-em-portugus)
9. [Referências](#9-referncias)
10. [Limitações explícitas](#10-limitaes-explcitas)
— Apêndices: [A. Tarefa → pacote/função em R](#apndice-a-tarefa--pacotefuno-em-r); [B. Glossário rápido](#apndice-b-glossrio-rpido); [C. Como rodar o script demonstrativo](#apndice-c-como-rodar-o-script-demonstrativo); [D. Comentários sobre exercícios](#apndice-d-comentrios-sobre-exercicios).

---

## 1. A ideia central: o que é um processo estocástico

### 1.1 Por que "processo", e não apenas "variável"?

Uma **variável aleatória** descreve a incerteza sobre *um* resultado: o valor de um dado lançado, a altura de uma pessoa sorteada, a demanda de amanhã. Quando a incerteza se desenrola **ao longo do tempo** — e cada instante pode produzir novos valores aleatórios —, precisamos de um objeto mais rico: um **processo estocástico**.

Pense nestes exemplos, todos ligados aos capítulos seguintes:

- O preço de uma ação ao longo dos meses: sobe, desce, sobe de novo… ninguém sabe exatamente o próximo passo, embora haja regularidades estatísticas.
- O número de clientes numa fila de banco entre meio-dia e uma da tarde.
- A concentração de um poluente num rio, medida em pontos e horários diferentes.
- O número de indivíduos infectados numa epidemia, dia após dia.
- A posição de uma molécula num fluido, agitada pelas colisões com outras moléculas (movimento browniano).

Em todos eles, existe um padrão estatístico subjacente — mas o detalhe individual é imprevisível. É exatamente esse tipo de fenômeno que a teoria dos processos estocásticos modela: "fenômenos que podem ser pensados como evoluindo no tempo de maneira aleatória" ([An Introduction to Stochastic Processes in Continuous Time](https://www.math.leidenuniv.nl/~spieksma/colleges/sp-master/LN.pdf)).

### 1.2 A definição formal — e o que ela diz na prática

> **Definição (processo estocástico).** Seja $(\Omega, \mathcal{F}, P)$ um espaço de probabilidade, $T$ um conjunto de índices (interpretado como "tempo") e $S$ um espaço de estados. Um **processo estocástico** com conjunto de parâmetros $T$ e espaço de estados $S$ é uma família de variáveis aleatórias
> $$X = (X_t)_{t \in T}, \quad X_t : \Omega \to S,$$
> todas definidas no mesmo espaço de probabilidade. Em linguagem coloquial: para cada instante $t$, $X_t$ é uma variável aleatória; juntas, essas variáveis compõem um sistema que evolui aleatoriamente.

Essa formulação aparece em livros-texto brasileiros, como nas notas do [Curso IPRJ/PUC-Rio](https://www.maxwell.vrac.puc-rio.br/15321/15321_4.PDF) e no capítulo introdutório do livro de [Introdução aos Processos Estocásticos (UFABC)](https://danielmiranda.prof.ufabc.edu.br/livros/processos/capitulo-01-processos-estocasticos.html), cuja notação seguimos aqui. A visão geral também está consolidada no verbete [Stochastic process](https://en.wikipedia.org/wiki/Stochastic_process): um objeto matemático geralmente definido como uma família de variáveis aleatórias indexadas, cujo índice costuma ser interpretado como tempo.

Traduzindo cada ingrediente:

- $\Omega$ (ômega): o conjunto de todos os resultados possíveis de toda a experiência ("o universo" das possibilidades).
- $\mathcal{F}$ (a álgebra de eventos): quais perguntas sobre o resultado fazemos e recebemos resposta "sim/não". Tecnicamente, conjuntos mensuráveis.
- $P$: a regra que atribui probabilidades aos eventos, com $P(\Omega)=1$.
- $T$: o conjunto de instantes observados. Pode ser discreto ($T = \{0, 1, 2, \dots\}$) ou contínuo ($T = [0, \infty)$).
- $S$: onde o processo vive. Pode ser finito (ex.: "saudável/doente"), enumerável (ex.: contagens $\{0, 1, 2, \dots\}$) ou contínuo (ex.: posições na reta).

Note que a palavra "processo" **não exige** que o índice seja temporal: um campo aleatório em uma imagem é indexado por coordenadas de pixels, e uma série de Fourier aleatória por frequências. Mas o uso esmagadoramente dominante é como modelo de evolução temporal.

### 1.3 Três exemplos construtivos (para fixar a ideia)

**(a) Sequência de moedas honestas.** Jogue uma moeda justa repetidamente. O espaço de resultados $\Omega$ pode ser pensado como o conjunto de sequências infinitas de caras ($H$) e coroas ($T$). Defina $X_n = 1$ se a $n$-ésima jogada for cara e $0$ caso contrário. Cada $X_n$ é uma cópia de um Bernoulli$(1/2)$, e juntos formam um processo estocástico em tempo discreto com espaço de estados $\{0, 1\}$. Este é o "átomo" de quase todos os processos: muitos modelos são sequências de decisões tipo cara/coroa com regras de transição mais elaboradas.

**(b) Posição de uma partícula na reta (passeio aleatório).** Uma partícula começa em $0$ e, a cada passo, anda $+1$ ou $-1$ com igual probabilidade, independentemente do passado. A posição no passo $n$ é $S_n = Z_1 + \cdots + Z_n$ com $Z_i \in \{-1, +1\}$. Aqui $T=\{0,1,2,\dots\}$ e $S=\mathbb{Z}$ (espaço de estados infinito enumerável). Esse processo será nosso primeiro grande personagem na seção 4.1.

**(c) Contagem de eventos até o instante $t$ (processo de contagem).** Considere chamadas telefônicas chegando a uma central. Para cada $t \geq 0$, defina $N(t)$ = "número de chamadas chegadas entre $0$ e $t$". Cada $N(t)$ é uma variável de contagem; o vetor inteiro $(N(t))_{t \geq 0}$ é um processo com tempo contínuo e espaço de estados $\{0, 1, 2, \dots\}$. Quando impomos axiomas razoáveis de independência e taxa constante, obtemos o **processo de Poisson** (seção 4.3).

### 1.4 Duas leituras que você precisa separar

Esta é a distinção conceitual mais importante da seção — e fonte de muita confusão:

1. **Instante fixo, resultado variável:** fixe $t$ e olhe $X_t$ como variável aleatória. Perguntamos: "qual a distribuição do estado no instante $t$?" Ex.: qual a chance de chover dia 10? (Ou: de haver 3 clientes na fila às 15h?)
2. **Resultado fixo, tempo variável:** fixe uma realização $\omega$ (o "sorteio" específico do universo) e olhe a função $t \mapsto X_t(\omega)$. Essa função é a **trajetória**, **caminho amostral** ou **realização** do processo. Ex.: a série histórica real de chuva do ano passado é uma trajetória.

Um gráfico de uma simulação mostra **uma trajetória**; o processo em si é o *mecanismo probabilístico* que gera todas as trajetórias possíveis. Você nunca vê o processo inteiro diretamente — só amostras dele. Como resumem as notas da UFABC: "um processo estocástico pode ser pensado como uma regra aleatória que escolhe uma trajetória; antes de sabermos qual ocorreu, temos várias trajetórias possíveis; depois que $\omega$ é fixado, uma delas fica determinada".

Quando quiser comunicar isso em uma frase: *a variável aleatória responde "onde podemos estar"; a trajetória conta "como foi desta vez"*.

### 1.5 De onde vêm os processos? A construção de Kolmogorov

Como se *constrói* um processo de fato? Uma abordagem moderna (devida a Andrey Kolmogorov) diz: especifique primeiro, para cada número finito de instantes $t_1, \dots, t_n$, a distribuição conjunta de $(X_{t_1}, \dots, X_{t_n})$ — as chamadas **distribuições de dimensão finita** —, e faça isso de forma **consistente**: a distribuição que você especifica para $n+1$ instantes deve ter como marginais as distribuições que você já tinha para $n$ instantes (a "condição de consistência"). Sob condições adequadas, existe então uma única lei de probabilidade que realiza todas essas distribuições simultaneamente — isto é, um processo com essas leis de dimensão finita. Esse é o conteúdo do [teorema de extensão de Kolmogorov](https://en.wikipedia.org/wiki/Kolmogorov_extension_theorem).

A intuição é poderosa: para modelar um processo, muitas vezes basta saber "que vetores de valores podem aparecer e com que probabilidade", e garantir que essas especificações não se contradigam quando vistas de forma marginal. Há sutilezas técnicas importantes (o teorema clássico garante a existência de um processo, mas não garante, por si só, regularidade das trajetórias — continuidade, por exemplo; isso exige controles adicionais, cf. [discussão no Risk Hub](https://riskhub.org/blogs/kolmogorovs-extension-perspective-on-stochastic-processes) e a literatura especializada). Na prática, usamos o teorema como garantia de que nossas especificações fazem sentido, e recorremos a construções explícitas (como limites de passeio aleatório) para obter trajetórias bem comportadas.

### 1.6 Classificação rápida (quadro de bolso)

| Tempo | Espaço de estados | Exemplo típico | Seção |
|---|---|---|---|
| Discreto | Finito | Cadeias de Markov em dois estados | 4.2 |
| Discreto | Enumerável | Passeio aleatório em $\mathbb{Z}$ | 4.1 |
| Discreto | Contínuo | Modelo AR(1) em cada passo | 4.7 |
| Contínuo | Enumerável | Processo de Poisson; cadeias de salto | 4.3, 4.2 |
| Contínuo | Contínuo | Movimento browniano; difusões | 4.4 |

Eixos adicionais úteis: **com/sem memória** (Markov versus processos com memória longa), **incrementos independentes** (Poisson, Wiener), **estacionariedade** (as estatísticas mudam com o tempo?), e **linearidade/não linearidade** das equações de evolução.

### Exercícios da seção 1

1. Dê três exemplos de processos estocásticos da sua área de interesse, indicando $T$ e $S$ em cada caso.
2. Explique, com suas palavras, a diferença entre "fixar $t$" e "fixar $\omega$". Onde moram os gráficos de simulação nessa dicotomia?
3. Uma urna tem bolas numeradas; sorteamos com reposição e anotamos a soma acumulada após $n$ sorteios. Descreva esse processo informalmente como família de variáveis aleatórias — qual é o espaço de estados?

*Soluções comentadas no Apêndice D.*


## 2. Ferramentas essenciais de probabilidade

Esta seção é uma revisão rápida e pragmática: as peças de probabilidade que reaparecem em todos os capítulos seguintes. Se você já domina esses tópicos, pode passar direto para a seção 3.

### 2.1 As distribuições que você vai encontrar em toda parte

| Distribuição | Notação | Onde aparece nos processos | Funções em R |
|---|---|---|---|
| Bernoulli | $X \sim \mathrm{Ber}(p)$ | cada passo de um passeio aleatório; chegada/não chegada num intervalo | `rbinom(1, 1, p)` |
| Binomial | $X \sim \mathrm{Bin}(n, p)$ | número de sucessos em $n$ ensaios; converte para Poisson quando $n$ grande, $p$ pequeno | `rbinom(n, n, p)` |
| Poisson | $X \sim \mathrm{Pois}(\lambda)$ | contagem de eventos raros num intervalo; número de chegadas por unidade de tempo | `rpois(n, lambda)` |
| Geométrica | $X \sim \mathrm{Geom}(p)$ | tempo até o primeiro sucesso (moedas até dar cara) | `rgeom(n, p)` |
| Exponencial | $X \sim \mathrm{Exp}(\lambda)$ | tempo entre chegadas no processo de Poisson; memória sem memória | `rexp(n, rate = lambda)` |
| Normal | $X \sim N(\mu, \sigma^2)$ | incrementos do movimento browniano; limites do TCL; ruído térmico | `rnorm(n, mean = mu, sd = sigma)` |
| Uniforme | $X \sim U(a, b)$ | incerteza "sem opinião"; simulação via inversão | `runif(n, min = a, max = b)` |

Duas conexões merecem ser guardadas:

1. **Poisson como limite do binomial.** Se $p_n \to 0$ e $n\,p_n \to \lambda > 0$, o número de sucessos converge em distribuição para $\mathrm{Pois}(\lambda)$. É a lei dos "eventos raros": por isso o Poisson modela terremotos, falhas e chegadas discretas com taxa média conhecida ([basic concepts of the Poisson process](https://www.probabilitycourse.com/chapter11/11_1_2_basic_concepts_of_the_poisson_process.php)).
2. **Exponencial e memória sem memória.** $X \sim \mathrm{Exp}(\lambda)$ satisfaz $P(X > x + a \mid X > a) = P(X > x)$ para $a, x \geq 0$: "o tempo restante não depende de quanto já passou". Essa propriedade caracteriza a exponencial (a menos da parametrização) e será decisiva na construção do Poisson (seção 4.3).

O acervo completo de densidades/funções de probabilidade do R segue o esquema `d*` (densidade), `p*` (acumulada), `q*` (quantil), `r*` (geração): `dbinom`, `pbinom`, `qbinom`, `rbinom`; `dpois`, `ppois`, …; `dnorm`, `pnorm`, `qnorm`, `rnorm`.

### 2.2 Esperança e variância (o essencial)

A **esperança** $E[X]$ é a média ponderada pelos valores possíveis: $\sum_x x\,P(X = x)$ (finitos ou somáveis). A **variância** $\mathrm{Var}(X) = E[(X - \mu)^2] = E[X^2] - \mu^2$ mede a dispersão típica ao quadrado, e seu expoente, o desvio-padrão $\sigma = \sqrt{\mathrm{Var}}$, mede dispersão na escala original.

Propriedades usadas o tempo todo:

- Linearidade: $E[aX + bY] = aE[X] + bE[Y]$; $\mathrm{Var}(aX) = a^2\,\mathrm{Var}(X)$.
- Soma de variáveis **independentes**: $\mathrm{Var}(X + Y) = \mathrm{Var}(X) + \mathrm{Var}(Y)$. Em geral há termos de covariância, e independência os elimina.
- Somando $n$ cópias i.i.d. de média $\mu$ e variância $\sigma^2$: $E[S_n] = n\mu$ e $\mathrm{Var}(S_n) = n\sigma^2$; logo o desvio-padrão de $S_n$ cresce como $\sqrt{n}$, enquanto a média empírica $S_n/n$ tem desvio-padrão $\sigma/\sqrt{n} \to 0$.

Esse último ponto é a semente dos dois teoremas-limite fundamentais: a concentração da média em torno de $\mu$ (lei dos grandes números) e a forma gaussiana da flutuação (teorema central do limite).

### 2.3 Independência e probabilidade condicional

Eventos $A$ e $B$ são **independentes** quando $P(A \cap B) = P(A)P(B)$, ou equivalentemente $P(A \mid B) = P(A)$: saber que $B$ ocorreu não altera a chance de $A$. Duas variáveis aleatórias são independentes quando qualquer evento gerado por $X$ é independente de qualquer evento gerado por $Y$.

Na modelagem estocástica, quase todas as hipóteses interessantes são afirmações de independência *condicional*: “dado o estado atual, o futuro não depende do passado” (Markov); “dados os instantes, os incrementos são independentes” (Poisson/Wiener); “dado o histórico até agora, a próxima jogada é justa” (martingala). Por isso vale entender bem a probabilidade condicional $P(A \mid B) = P(A \cap B)/P(B)$ e sua versão para variáveis, $E[X \mid Y]$, tratada em 2.6.

### 2.4 Lei dos Grandes Números (LGN)

Em linguagem simples: se repetimos um experimento aleatório $n$ vezes de forma independente, a **média empírica** dos resultados se aproxima da esperança teórica conforme $n$ cresce. Para moedas honestas, a frequência de caras tende a $1/2$; para dados honestos, a média dos lançamentos tende a $3{,}5$.

Uma versão simulável em poucos minutos de R (veja também §5):

```{r}
set.seed(42)
n <- 200                 # tamanho de cada amostra
nsim <- 1000             # número de repetições do experimento
medias <- replicate(nsim, mean(sample(1:6, size = n, replace = TRUE)))
mean(medias)             # deve estar perto de 3.5
sd(medias)               # ≈ sqrt(35/12)/sqrt(n) ≈ 0.17
```

Você deve observar que a distribuição das médias se concentra cada vez mais em torno de $3{,}5$ conforme $n$ cresce — a variabilidade da média cai como $1/n$, um padrão tão robusto que fundamenta praticamente toda a inferência estatística ([Law of large numbers and the Central Limit Theorem — notas do MIT](https://math.mit.edu/~dav/05.dir/class6-prep.pdf)).

### 2.5 Teorema Central do Limite (TCL)

O TCL responde à pergunta natural: qual é a *forma* da flutuação da média em torno de $\mu$ para amostras grandes? Resposta: aproximadamente **Normal**. Mais precisamente, para $S_n = X_1 + \cdots + X_n$ com $X_i$ i.i.d., média $\mu$ e variância $\sigma^2 < \infty$:

$$\frac{S_n - n\mu}{\sigma \sqrt{n}} \;\xrightarrow{d}\; N(0, 1).$$

Intuicionamento: mesmo que os dados venham de uma distribuição muito assimétrica (ex.: uniforme discreta de um dado), a média padronizada vira gaussiana quando muitos termos são somados. Isso vale ouro porque (i) permite intervalos de confiança e testes sem conhecer a distribuição exata de cada fenômeno; (ii) explica por que a Normal aparece em tantas modelagens — incluindo como incremento do movimento browniano (seção 4.4). Uma boa intuição dinâmica está nas [notas do MIT sobre LGN e TCL](https://math.mit.edu/~dav/05.dir/class6-prep.pdf); o enunciado preciso e condições de validade estão no verbete [Central limit theorem](https://en.wikipedia.org/wiki/Central_limit_theorem).

Simulação que mostra o fenômeno em ação:

```{r}
set.seed(123)
# Amostragem da distribuição de um único dado (uniforme em {1,...,6})
x1 <- sample(1:6, 20000, replace = TRUE)
# Médias de amostras de tamanhos crescentes: a curva fica Normal
par(mfrow = c(1, 3))
for (n in c(2, 10, 100)) {
  mm <- replicate(5000, mean(sample(x1, n, replace = TRUE)))
  hist(mm, breaks = 40, main = paste("n =", n), col = "gray80", border = "white")
  abline(v = mean(x1), col = "red", lwd = 2)
}
```

Compare a figura com a intuição: para $n=2$, a distribuição das médias ainda carrega o formato do dado; para $n=100$, é indistinguível de uma Normal.

### 2.6 Esperança condicional: a ponte para as martingalas

Considere duas variáveis $X$ e $Y$. A **esperança condicional** de $X$ dado $Y$ é a função (v.a.) $\omega \mapsto E[X \mid Y](\omega)$ que, para cada valor observado de $Y$, dá a média de $X$ sob a distribuição de $X$ *atualizada* por aquela informação. Propriedades úteis e válidas em generalidade:

- $E[X \mid Y]$ é uma **variável aleatória**, função de $Y$; $E[E[X \mid Y]] = E[X]$ (propriedade da torre: suavizar duas vezes é o mesmo que suavizar uma).
- Se $X$ é independente de $Y$, então $E[X \mid Y] = E[X]$ (informação irrelevante).
- Desigualdade de Jensen: para $f$ convexa, $E[f(X) \mid Y] \geq f(E[X \mid Y])$ (bom saber: aparecerá em argumentos de convergência de martingalas).

Exemplo pedagógico (e precursor da seção 3): no passeio aleatório simétrico, dado que hoje a posição é $k$, o valor esperado da posição amanhã é exatamente $k$ — porque o próximo passo é $+1$ ou $-1$ com igual probabilidade: $E[S_{n+1} \mid S_n = k] = k$. Esse “não ganho nem perco na expectativa, com toda informação passada” é a marca registrada dos **processos de jogo justo**, isto é, martingalas.

### Exercícios da seção 2

1. Simule $10^4$ jogadas de um dado honesto e compare a frequência de cada face com $1/6$. O que muda se você usar $10^6$ jogadas?
2. Explique por que, numa urna com reposição, o vetor $(\text{total após } n \text{ sorteios})_{n \geq 0}$ tem incrementos independentes.
3. Seja $X \sim \mathrm{Exp}(2)$ (taxa). Calcule analiticamente $P(X > 1)$ e simule para conferir. Qual a memória sem memória dizendo $P(X > 2 \mid X > 1)$?

*Soluções comentadas no Apêndice D.*


## 3. Informação no tempo: filtrações, tempos de parada e martingalas

Esta seção formaliza uma ideia que já está implícita em tudo o que fizemos até aqui: "o que sabemos quando decidimos". Ela culmina no conceito de **martingala**, a classe mais elegante de processos e a chave técnica para análise de limites, testes sequenciais e finanças quantitativas.

### 3.1 A informação cresce — a noção de filtração

Ao observar um processo ao longo do tempo, sua base de conhecimento só aumenta: depois de 10 jogadas da moeda, você sabe tudo sobre as 10 primeiras; depois de 20, também sobre as 10 primeiras (nada foi apagado). Essa estrutura axiomática é capturada pela ideia de **filtração**.

Dado um processo $X = (X_t)_{t \in T}$ definindo seu espaço amostral, o objeto matemático que representa "as perguntas respondíveis até o instante $t$" é uma álgebra $\sigma$ de eventos, notadamente

$$\mathcal{F}_t = \sigma\{X_s : s \leq t\},$$

a menor álgebra que torna todas as variáveis $X_s$ com $s \leq t$ mensuráveis. A família $(\mathcal{F}_t)_{t \in T}$ tem duas propriedades:

- **crescente**: se $s \leq t$, então $\mathcal{F}_s \subseteq \mathcal{F}_t$ (informação não se perde);
- **contida** na álgebra final $\mathcal{F}$.

Quando uma segunda fonte de informação for maior — por exemplo, além das cotações $X_t$ você observa também notícias $Y_t$ —, usamos famílias $\sigma$ maiores, ainda crescendo com $t$. Em resumo, seguindo as notas de [filtrations and stopping times](https://www.randomservices.org/random/prob/Stop.html):

> Uma **filtração** é simplesmente uma família crescente de subálgebras $\sigma$, indexada pelo tempo; pensamos em $\mathcal{F}_t$ como a coleção de eventos cujo status ("ocorreu ou não") pode ser decidido observando o processo até $t$.

Um processo $Y = (Y_t)$ é dito **adaptado** à filtração $(\mathcal{F}_t)$ quando, para cada $t$, a variável $Y_t$ é $\mathcal{F}_t$-mensurável — isto é, "o valor de $Y_t$ é compatível com a informação disponível até $t$": dá para saber se $Y_t$ caiu num conjunto medindo apenas o passado até $t$. Se isso vale para toda trajetória razoável, o processo é adaptado; quase todo processo construído "causalmente" é.

### 3.2 Tempos de parada: parar sem olhar o futuro

Uma regra de decisão sobre *quando* agir (parar de apostar, encerrar um experimento, vender um ativo) deve ser tomada apenas com base no passado e no presente: nada de usar o futuro. Matematicamente:

> **Definição (tempo de parada).** Um tempo aleatório $T$ tomando valores em $T \cup \{\infty\}$ é um **tempo de parada** (relativo à filtração $(\mathcal{F}_t)$) quando, para todo $t$, o evento $\{T \leq t\}$ pertence a $\mathcal{F}_t$ — ou seja, é possível decidir "já parei até $t$?" olhando apenas a informação acumulada até $t$.

Exemplos importantes:

- **Tempo de primeira visita.** Para $a > 0$, $\tau_a = \inf\{t : X_t \geq a\}$ é um tempo de parada (nos casos discretos que usaremos, basta poder decidir "já atingi $a$?" olhando apenas o histórico até $t$; para difusões em tempo contínuo vale sob hipóteses de continuidade adequadas). Decidir "pare quando eu atingir lucro $a$".
- **Vencimento de um prazo.** Fixar datas determinísticas ($T = 5$ dias) é trivialmente um tempo de parada — decisões pré-comprometidas são permitidas.
- **Primeira cara numa sequência de moedas.** $T$ = número da jogada onde sai cara pela primeira vez. É um tempo de parada com distribuição geométrica.
- **Contraexemplo clássico.** O tempo do *último* retorno à origem antes de um instante fixo $N$ (ex.: "qual foi o último dia em que a partícula esteve em $0$, olhando toda a trajetória até hoje?") **não** é tempo de parada: para saber se aquele foi o último retorno, é necessário conhecer o futuro dentro do intervalo ([notes de aula CMU, "Random Times and Their Properties"](https://www.stat.cmu.edu/~cshalizi/754/notes/lecture-06.pdf)).

A distinção entre tempos de parada e tempos "olhadores de futuro" é sutil mas crítica: ela separa os argumentos probatórios válidos dos paradoxos clássicos.

### 3.3 Martingalas: o jogo justo idealizado

Vamos agora ao protagonista desta seção. Intuitivamente, uma martingala é um modelo de **jogo justo**: dado todo o histórico, a expectativa futura do capital é exatamente o capital atual — nem vantagem, nem desvantagem sistemática. Historicamente, a ideia veio de estratégias de apostas que anulam perdas (o nome vem do inglês *martingale*, uma família de estratégias de dobrar a aposta), mas a noção moderna é bem mais ampla e bela.

> **Definição (martingala).** Seja $(\Omega, \mathcal{F}, (\mathcal{F}_t)_{t \in T}, P)$ um espaço probabilístico filtrado. Um processo real $M = (M_t)_{t \in T}$ é uma **martingala** se:
> 1. $M_t$ é integrável ($E[|M_t|] < \infty$) para todo $t$;
> 2. $M$ é adaptado a $(\mathcal{F}_t)$;
> 3. para todos $s \leq t$: $E[M_t \mid \mathcal{F}_s] = M_s$ (quase certamente).
>
> Equivalentemente, para incrementos: $E[M_{t+1} - M_t \mid \mathcal{F}_t] = 0$ — o próximo passo tem expectativa zero, sabendo tudo até agora.

Observações úteis:

- **Submartingala/supermartingala:** se a condição 3 vira desigualdade $E[M_t \mid \mathcal{F}_s] \geq M_s$ (resp. $\leq$), diz-se sub (resp. super). Sub + super = martingala; o sinal indica tendência direcional do processo ("jogo com vantagem" ou "desvantagem").
- **Localização (local martingales):** algumas propriedades interessantes requerem condições que valem apenas "até tempos de parada crescentes"; define-se então a classe dos *martingalas locais*, importante em cálculo estocástico (por exemplo, integrais estocásticas são localmente martingalas).
- Toda martingala satisfaz a chamada **propriedade de média**: o valor atual é a média ponderada de qualquer distribuição futura consistente com a informação presente.

Exemplos centrais que aparecerão adiante:

1. **Passeio aleatório simétrico** é martingala: $E[S_{n+1} \mid S_n] = S_n$ (cada passo $Z_i$ tem média zero independente do passado).
2. **Passeio aleatório enviesado NÃO é** martingala: $E[Z_i] = 2p - 1 \neq 0$ gera deriva sistemática ($\pm 1$ com probabilidades diferentes). Você pode compensar a deriva subtraindo $n(2p-1)$, recuperando um processo martingálico.
3. **Movimento browniano** $W_t$ é martingala quadrado-integrável com $E[W_t^2] = t$ (seção 4.4).
4. **Processos de Poisson centrados:** $N(t) - \lambda t$ é martingala (subtrair a taxa média remove a tendência de crescimento).
5. **Wealth de um apostador:** sob odds justas (esperança nula de lucro por aposta), o capital ao longo do tempo é uma martingala. Se o cassino tem vantagem (jogo desonesto), é supermartingala — o capital tende a cair em expectativa. Este ponto de vista é explorado com profundidade em [Martingales and the Optional Stopping Theorem](https://www.jeremykun.com/2014/03/03/martingales-and-the-optional-stopping-theorem/).

### 3.4 O teorema da parada opcional (Doob)

A pergunta prática: em um jogo justo, se você decide parar num momento adaptado (digamos, ao bater a meta), sua expectativa final ainda é igual à inicial? Resposta em condições adequadas: **sim**. Esse é o conteúdo do **teorema da parada opcional de Doob**:

> **Teorema (Doob, forma simplificada).** Sejam $M$ uma martingala e $T$ um tempo de parada. Se $M_T$ é integrável e uma das hipóteses técnicas adicionais vale (por exemplo, $T$ é limitado, ou a família $(M_t)$ é uniformemente integrável, ou $M_T \geq |M_0|$... variantes existem), então $E[M_T] = E[M_0]$.

Em linguagem de apostas: nenhuma estratégia de parada adaptada transforma um jogo justo em máquina de lucro esperado; quem tenta geralmente ignora alguma condição técnica e acaba surpreendido (paradoxos famosos envolvem estratégias agressivas de dobrar a aposta (tipo "doubling strategy", que violam integrabilidade uniforme). O outro lado igualmente importante: quando o jogo tem desvio constante (super/submartingala), a desigualdade se propaga às paradas opcionais — é o que sustenta, por exemplo, que casa de apostas lucra a longo prazo sob certas condições e que certos jogos de azar levam à ruína com probabilidade alta. Veja discussões acessíveis em [Jeremy Kun's primer on martingales](https://www.jeremykun.com/2014/03/03/martingales-and-the-optional-stopping-theorem/) e nas [notes of Marc Yor / Samuel Drapeau](https://www.samuel-drapeau.info/downloads/WS2015/StochProcess/Stochastic_Processes_chap02_1126.pdf).

Aplicações do OST que veremos depois:

- **Ruína do jogador:** calcular probabilidades de absorção em barreiras usando martingalas apropriadas (seção 4.1).
- **Testes sequenciais:** planejamento de ensaios com amostragem até decisão (revisões Bayesianas usam essa ideia continuamente).
- **Finanças:** precificação de derivativos americanos e estratégias de hedging ótimas em modelos de mercado com parada ótima.

### 3.5 Onde isso reaparece neste guia

| Conceito | Reaparece em | Papel |
|---|---|---|
| Filtração | 4.x, 6.x | definir adaptação, causalidade, informação |
| Tempo de parada | 4.1, 4.4, 6.A, 7 | parar em metas, primeiro retorno, exercício de opção |
| Martingala | 4.4, 6.A, 8 | passeio simétrico, BM, precificação justa |
| Parada opcional | 4.1, 6.A | ruína do jogador, opções americanas, validação de simulações |

### Exercícios da seção 3

1. Mostre diretamente que $S_n = \sum Z_i$ com $Z_i$ i.i.d. e $E[Z_i] = 0$ é martingala. Por que a mesma prova falha para o passeio enviesado?
2. Verifique se $N(t) - \lambda t$ (processo de Poisson centrado) é martingala — use a definição e as propriedades de incrementos independentes do processo de Poisson (você pode conferir na seção 4.3).
3. Explique (sem tecnicismos) por que a última passagem pela origem não é um tempo de parada.
4. No problema da moeda honesta até sair cara (geométrico), escreva a recursão para $E[T]$ e resolva — compare com a intuição de que $E[T] = 1/p$ (veja [ABRACADABRA problem via martingales](https://www.jeremykun.com/2014/03/03/martingales-and-the-optional-stopping-theorem/)).

*Soluções comentadas no Apêndice D.*


## 4. As grandes famílias de processos

Esta seção é o coração do guia. Para cada classe importante de processo estocástico, seguimos o mesmo roteiro: (1) intuição e definição; (2) propriedades-chave; (3) como simular em R. A ideia é que você reconheça o "molde" de cada modelo — e saiba qual usar quando.

#> **Figura de referência (simulações):** ![](doc/figs/exemplos_processos.png)

## 4.1 Passeio aleatório simples (o movimento browniano discreto)

**Intuição.** Imagine uma pessoa embriagada caminhando numa calçada retilínea, dando passos de tamanho fixo para a frente ou para trás, com igual probabilidade, sem memória direcional. A posição após $n$ passos é uma soma de $\pm 1$ independentes. Esse é o passeio aleatório simples — o processo estocástico mais fundamental que existe, e o ponto de partida de quase toda a teoria (a partir dele se constroem aproximações para difusões, séries temporais, modelos de preços etc.).

**Definição formal.** Sejam $Z_1, Z_2, \dots$ variáveis i.i.d. com $P(Z_i = +1) = p$ e $P(Z_i = -1) = q = 1 - p$. O **passeio aleatório** é o processo $S_n = Z_1 + \cdots + Z_n$ (com $S_0 = 0$). O caso $p = q = 1/2$ é o *passeio simétrico*; $p \neq 1/2$ dá um *passeio enviesado* (com deriva).

Cálculos elementares mas onipresentes (guarde-os):

$$E[S_n] = n\,(p - q) = n(2p - 1), \qquad \mathrm{Var}(S_n) = 4n\,p\,q,$$

pois $E[Z_i] = p - q$ e $\mathrm{Var}(Z_i) = 1 - (p-q)^2 = 4pq$. Em particular, no caso simétrico: $E[S_n] = 0$ e $\mathrm{Var}(S_n) = n$ — o desvio-padrão cresce como $\sqrt{n}$, sinal de que a partícula tipicamente se encontra a uma distância da ordem de $\sqrt{n}$ da origem após $n$ passos (não $n$! os passos para lá e para cá se cancelam em parte).

**Recorrência versus transiência.** Duas perguntas profundas: a partícula volta à origem infinitas vezes (recorrente) ou apenas finitamente (transiente)? Aqui surge um dos resultados mais notáveis da teoria, o **teorema de Pólya**:

- O passeio simétrico em dimensões $d = 1$ e $d = 2$ é **recorrente**: retorna à origem com probabilidade $1$ (e passa por ela infinitas vezes).
- Em dimensões $d \geq 3$, é **transiente**: há probabilidade positiva de nunca mais voltar à origem — em três dimensões, a partícula consegue "se perder".

O resultado clássico foi formulado por George Pólya em 1919–1921 ([notas introdutórias](https://www.cmat.edu.uy/~lessa/resource/randomwalknotes.pdf), [resumo didático](https://math.dartmouth.edu/~pw/math100w13/mare.pdf)); discussões intuitivas podem ser vistas em [por que em 1D/2D a recusa é certa](https://math.stackexchange.com/questions/3233552/intuitive-reason-why-a-simple-symmetric-random-walk-is-recurrent-on-bbb-z2-a) e em [Por que passeios aleatórios se perdem em 3D](https://www.youtube.com/watch?v=byvEzyFgv44). Uma pista heurística: em dimensão baixa, os "caminhos que voltam" são numerosos o bastante para dominar; em dimensão alta, o volume de escape cresce e a soma das probabilidades de retorno converge — a integral que decide a recorrência diverge só em dimensões baixas.

```{r}
# --- Simulação do passeio aleatório simples em 1D -------------------------
set.seed(42)
n <- 300                                   # número de passos
passos <- sample(c(-1, 1), n, replace = TRUE)
caminho <- c(0, cumsum(passos))            # posição S_0, S_1, ..., S_n

plot(caminho, type = "l", col = "steelblue", lwd = 1.8,
     xlab = "passo", ylab = "posição", main = "Passeio aleatório simétrico em 1D (300 passos)")
abline(h = 0, col = "gray60")
lines(caminho, col = "steelblue")
points(seq_along(caminho), caminho, pch = 20, cex = 0.7, col = "navy")
```

Você deve observar: (i) a trajetória oscila de forma irregular, sem tendência; (ii) as faixas típicas ficam dentro de $\pm C\sqrt{n}$; (iii) quanto mais longa a janela observada, mais retornos à vizinhança da origem aparecem (em dimensão 1 e 2 isso persiste indefinidamente — enquanto em 3D você verá fugas definitivas).

Para experimentar dimensões maiores (simulação ilustrativa, muitos passos):

```{r}
# --- Retorno à origem em dimensões 1, 2 e 3 -------------------------------
set.seed(7)
sim_rw_dim <- function(nsteps, dim){
  steps <- matrix(sample(c(-1, 1), nsteps * dim, replace = TRUE), nrow = nsteps, ncol = dim)
  pos   <- rbind(rep(0, dim), apply(steps, 2, cumsum))     # inclui o tempo 0
  return(pos)                                              # cada linha = posição num instante
}
conta_retornos <- function(pos, dim){
  zero <- rowSums(pos == 0) == dim                         # todos os componentes nulos?
  sum(zero[ -1 ])                                          # ignorar o instante inicial
}

nsim <- 500; nsteps <- 2000
resultados <- lapply(1:3, function(d){
  replicate(nsim, conta_retornos(sim_rw_dim(nsteps, d), d))
})
sapply(resultados, function(x) mean(x > 0))   # fração de simulações com algum retorno
# Esperado qualitativamente: ~1.0 (d=1), ~1.0 (d=2), < 1 (d=3)
```

*(Aviso: rodar isso consome alguns segundos; reduza `nsim` para começar.)*

**Ruína do jogador (gambler's ruin).** Um dos problemas históricos mais importantes: dois apostadores com patrimônios iniciais $i$ e $N - i$ apostam repetidamente $1$ unidade por rodada; o jogo termina quando um deles ficar sem dinheiro. Se o jogo for **justo** ($p = q = 1/2$), a probabilidade de o primeiro jogador atingir $N$ antes de ruir é simplesmente

$$P_i(\text{vitória}) = \frac{i}{N},$$

enquanto no caso **desonesto** ($p \neq q$):

$$P_i(\text{vitória}) = \frac{(q/p)^i - (q/p)^N}{1 - (q/p)^N}.$$

Simbólico: basta o capital inicial dividido pelo capital total no jogo justo. Se o jogo tem ligeira desvantagem, o expoente faz mágica ruim: mesmo começando com metade do capital, a probabilidade de vitória cai drasticamente com $N$ — e a mesa geralmente vence. É a versão matemática do ditado "a casa sempre ganha". (Ver [Gambler's ruin — Wikipédia](https://en.wikipedia.org/wiki/Gambler%27s_ruin), [nota histórica](http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-GR.pdf) e [derivação via funções geradoras, Grinstead & Snell](https://stats.libretexts.org/Bookshelves/Probability_Theory/Introductory_Probability_(Grinstead_and_Snell)/12%3A_Random_Walks/12.02%3A_Gambler%27s_Ruin).)

```{r}
# --- Ruína do jogador: fórmula vs. simulação -----------------------------
prob_teorica <- function(i, N, p = 0.5){
  if (abs(p - 0.5) < 1e-12) return(i / N)
  (((1 - p) / p)^i - ((1 - p) / p)^N) / (1 - ((1 - p) / p)^N)
}
ruina_sim <- function(i, N, p = 0.5, rodadas = 5000, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  vitorias <- 0
  for (r in 1:rodadas){
    saldo <- i
    while (saldo > 0 && saldo < N){
      saldo <- saldo + sample(c(-1, 1), size = 1, prob = c(q = 1 - p, p = p))
    }
    vitorias <- vitorias + (saldo == N)
  }
  c(simulacao = vitorias / rodadas, teorio = prob_teorica(i, N, p))
}
ruina_sim(i = 10, N = 50, p = 0.5, seed = 1)   # ambos ≈ 0.20
```

**Passeio aleatório como ponte para outros processos.** Quando passamos ao limite de passos muito pequenos e frequentes (escala adequada), o passeio aleatório simétrico converge para o **movimento browniano** — a versão contínua que veremos na seção 4.4. Essa conexão (teorema de Donsker–Skorokhod, informalmente) é o motivo pelo qual os dois objetos compartilham tantas propriedades: martingala, ausência de memória, autossemelhança.

**Aplicações típicas:** modelagem de difusão molecular (Einstein, 1905), estratégia de busca em redes, genética populacional (deriva genética), preços de ativos sob a hipótese de mercado eficiente, filas e tráfego.

### 4.2 Cadeias de Markov (memória que cabe numa matriz)

**Intuição.** Em muitos sistemas, o futuro depende do presente, mas não do resto do passado dado o presente: amanhã eu decido baseado em como estou hoje, não em como estava anteontem. Sistemas com essa propriedade são modelados por **cadeias de Markov**, cuja evolução é governada por uma matriz de transição.

**Definição (tempo discreto).** Um processo $(X_n)_{n \geq 0}$ com valores em $S$ (finito ou enumerável) é uma **cadeia de Markov** quando existe uma matriz $P = (p_{ij})$ de probabilidades de transição tal que

$$P(X_{n+1} = j \mid X_0 = i_0, \ldots, X_n = i) = p_{ij} \quad \text{para todo } j,$$

ou seja, a distribuição do próximo estado depende apenas do estado atual. Assumimos $\sum_j p_{ij} = 1$ para cada $i$ ("linha estocástica"). A identidade de **Chapman–Kolmogorov** expressa a consistência: $P^n = P \cdot P^{n-1}$, e seus elementos dão as probabilidades de ir de $i$ a $j$ em exatamente $n$ passos:

$$P^n(X_n = j \mid X_0 = i) = (P^n)_{ij}.$$

Exemplos clássicos: o clima de amanhã (seco/chuvoso) depende apenas do clima de hoje; a composição de uma espécie em equilíbrio mutacional (alelos com taxas de conversão fixas); o valor de uma carteira que muda segundo um ambiente de mercado em dois estados (alta/baixa); o famoso "problema dos elevadores"; o PageRank dos mecanismos de busca (com modificação técnica de teleporte).

**Classificação e comportamento assintótico.** Finitas cadeias conectadas exibem comportamento universal: depois de muito tempo, a distribuição de estados converge para uma distribuição estacionária única — independente da condição inicial. As condições técnicas envolvem **irredutibilidade** (qualquer estado alcançável de qualquer outro) e **aperiodicidade**; se houver absorção (estados que "aprisionam"), a convergência ocorre para uma combinação ponderada destes. A ferramenta central é a **medida estacionária** $\pi$: vetor de probabilidades que satisfaz

$$\pi P = \pi \quad\Longleftrightarrow\quad \sum_i \pi_i p_{ij} = \pi_j \;\; \text{para todo } j,$$

normalizada por $\sum_i \pi_i = 1$. Interpretando $\pi_i$ como proporção de longo prazo do tempo gasto no estado $i$ (ergodicidade da cadeia finita). A existência/unicidade sob irredutibilidade e aperiódico é o teorema ergódico de Perron–Frobenius para matrizes estocásticas ([tutorial espectral para cientistas de dados](https://arxiv.org/pdf/2207.02296)); a resolução prática costuma vir do autovetor associado ao autovalor $1$ ([estações estacionárias e limites](https://www.probabilitycourse.com/chapter11/11_3_2_stationary_and_limiting_distributions.php)).

**Tempo contínuo.** Uma cadeia de Markov a tempo contínuo (CTMC) tem matriz de intensidades de transição $Q = (q_{ij})$ (as taxas instantâneas), com equações forward/backward análogas; exemplos incluem epidemias SIR (cada indivíduo infectado transmite a taxa $\beta$ e recupera a taxa $\gamma$), nascimentos e mortes, e o processo de Poisson homogêneo como CTMC com apenas transições $i \to i+1$ de taxa $\lambda$ constante. Exemplos e formulações acessíveis estão nas ementas brasileiras (p.ex., [ME501 — Unicamp](https://www.ime.unicamp.br/graduacao/disciplina/me501-processos-estocasticos)) e nos textos clássicos.

```{r}
# --- Cadeia de Markov de 3 estados: distribuição estacionária ------------
P <- rbind(
  c(0.7, 0.2, 0.1),   # de A para (A, B, C)
  c(0.2, 0.6, 0.2),   # de B
  c(0.1, 0.3, 0.6)    # de C
)
rownames(P) <- colnames(P) <- c("A", "B", "C")

distrib_estacionaria <- function(Pmat){
  eig <- eigen(t(Pmat))              # autovetores-linha da transposta
  k   <- which.min(abs(eig$values - 1))   # autovalor mais próximo de 1
  pi  <- Re(eig$vectors[, k])
  names(pi) <- rownames(Pmat)
  pi / sum(pi)                       # normalizar
}
pi_hat <- distrib_estacionaria(P)
pi_hat           # ≈ (0.345, 0.379, 0.276)
pi_hat %*% P     # conferir: deve dar o próprio pi_hat
```

```{r}
# --- Simulação de trajetórias e frequências empíricas --------------------
sim_cadeia <- function(Pmat, n, semente = NULL, ini = 1){
  if (!is.null(semente)) set.seed(semente)
  K <- nrow(Pmat)
  estados <- integer(n)
  atual <- ini
  for (t in 1:n){
    estados[t] <- atual
    atual <- sample(K, size = 1, prob = Pmat[atual, ])
  }
  factor(estados, levels = seq_len(K), labels = rownames(Pmat))
}

set.seed(11)
traj <- sim_cadeia(P, n = 20000, ini = 1)
table(traj) / length(traj)    # frequência empírica ≈ pi_hat
head(traj, 20)                # primeiros 20 estados da trajetória
```

Você deve ver as frequências empíricas coincidindo (até erro Monte Carlo) com $\pi$ — e note que a cadeia "esquece" a condição inicial. Essa é a propriedade que torna as cadeias de Markov tão úteis em computação: elas convergem para sua distribuição estacionária, e muitas tarefas (amostragem de Boltzmann, ranking de páginas, análise de filas) usam exatamente essa convergência.

### Exercícios da seção 4

1. Verifique analiticamente que, no passeio simétrico 1D, $E(S_{2n}) = 0$ e que o número esperado de retornos à origem cresce com $n$ (dica: use a convolução das leis e a Stirling para estimativas assintóticas).
2. Estude o problema da ruína do jogador com $p = 0{,}55$ (favorável ao jogador!) e $i = 50$, $N = 100$. Compare o resultado exato com a intuição "jogo favorável" e explique matematicamente.
3. Numa cadeia com três estados, faça a matriz $P$ ter estados 1 e 2 recorrentes (duas classes fechadas disjuntas). Descreva o que acontece com a distribuição limite dependendo do estado inicial. Por que a unicidade de $\pi$ falha aqui?
4. Implemente uma versão vetorial da simulação da cadeia (`replicate` + vetorização) e compare a velocidade com a versão em loop — discuta o trade-off legibilidade × desempenho.

*Soluções comentadas no Apêndice D.*


### 4.3 Processo de Poisson (chegadas aleatórias com taxa constante)

**Intuição.** Suponha que eventos "raros e independentes" ocorram a uma taxa média conhecida: chamadas numa central, falhas num equipamento, fótons chegando a um sensor, desintoxicações em um hospital. O processo de Poisson formaliza a ideia de "aleatoriedade pura com frequência média fixa": os tempos entre eventos são aleatórios (exponenciais), mas o número esperado por unidade de tempo é estável.

**Definição (via incrementos).** Um processo de contagem $N = (N(t))_{t \geq 0}$ é um **processo de Poisson homogêneo de taxa $\lambda > 0$** quando:

1. $N(0) = 0$;
2. $N$ tem **incrementos independentes**: para intervalos disjuntos, os números de eventos são independentes;
3. $N$ tem **incrementos estacionários**: o número de eventos num intervalo depende apenas do comprimento do intervalo;
4. Para cada $t$, $N(t) \sim \mathrm{Pois}(\lambda t)$ — equivalently, $P(N(t) = k) = e^{-\lambda t}\frac{(\lambda t)^k}{k!}$.

Condições análogas modelam eventos no espaço (processo pontual de Poisson), onde "tempo" vira "área/volume" ([Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process)).

**Definição equivalente (interchegadas).** O processo também pode ser construído exigindo que os tempos entre chegadas sucessivas $T_k - T_{k-1}$ sejam i.i.d. $\sim \mathrm{Exp}(\lambda)$. A conexão com a memória sem memória da exponencial é o coração do modelo: "o processo não envelhece"; a probabilidade de esperar mais tempo não depende de quanto já se esperou ([basic concepts of the Poisson process](https://www.probabilitycourse.com/chapter11/11_1_2_basic_concepts_of_the_poisson_process.php); para um tratamento clássico, [MIT 6.262, capítulo 2](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/3a19ce0e02d0008877351bfa24f3716a_MIT6_262S11_chap02.pdf)).

Propriedades úteis:

- $E[N(t)] = \lambda t$ e $\mathrm{Var}(N(t)) = \lambda t$ (igualdade média-variância é assinatura do Poisson).
- **Superposição:** se somarmos dois processos de Poisson independentes de taxas $\lambda_1$ e $\lambda_2$, obtemos novamente Poisson de taxa $\lambda_1 + \lambda_2$ (atende muito bem a "fontes" aditivas de eventos).
- **Enfraquecimento ("thinning"):** se cada evento é mantido com probabilidade $p$ (independentemente de tudo), o processo dos eventos mantidos é Poisson de taxa $\lambda p$.
- **Não-homogêneo:** permitindo taxa $\lambda(t)$ dependente do tempo, obtém-se o processo de Poisson não-homogêneo, com $E[N(t)] = \int_0^t \lambda(s)\,ds$ — usado para modelos com sazonalidade ou tendências simples (por exemplo, ocorrências de raios durante uma tempestade).

```{r}
# --- Simulação do processo de Poisson ------------------------------------
set.seed(9)
lambda <- 3          # eventos por hora
Tmax   <- 4          # janela de observação
chegadas <- cumsum(rexp(200, rate = lambda))   # interchegadas exponenciais
chegadas <- chegadas[chegadas <= Tmax]         # manter dentro da janela

# Trajetória N(t) (função degrau)
tgrid <- seq(0, Tmax, by = 0.01)
Nt <- sapply(tgrid, function(t) sum(chegadas <= t))
plot(tgrid, Nt, type = "s", col = "darkgreen", lwd = 2,
     xlab = "tempo", ylab = "N(t)", main = "Processo de Poisson — trajetória simulada")
abline(v = chegadas, col = "gray70", lty = 3)

# Distribuição de N(t) num instante fixo vs. Poisson teórico
tt <- 2
nsim <- 3000
contagem <- replicate(nsim, {
  ch <- cumsum(rexp(300, rate = lambda))
  sum(ch <= tt)
})
hist(contagem, breaks = 0:14, freq = FALSE, col = "gray85", border = "white",
     main = paste("Distribuição empírica de N(", tt, ") vs. Poisson(", lambda*tt, ")", sep = ""))
points(0:14, dpois(0:14, lambda = lambda * tt), col = "red", pch = 19)
c(media_empirica = mean(contagem), var_empirica = var(contagem), teoria = lambda * tt)
```

Você deve observar: a média e a variância empíricas próximas de $\lambda t$ (linhas vermelhas batendo com as barras), e a trajetória com saltos unitários — característica típica de processos pontuais discretos.

**Aplicações típicas:** filas (input M/M/1 da seção 6.C), confiabilidade (falhas raras), epidemiologia (casos esporádicos), cosmologia (fótons), seguros (sinistros raros), e como "bloco de construção" de processos mais complexos (renovação, composto de Poisson).

### 4.4 Movimento browniano / processo de Wiener (o caos contínuo)

**História em três atos.** (1) Em 1827, Robert Brown observou grãos de pólen dançando na água ao microscópio — movimento errático aparentemente eterno. (2) Em 1905, Albert Einstein explicou-o como o efeito de colisões moleculares incessantes, dando previsões quantitativas verificáveis. (3) Na década de 1920, Norbert Wiener deu a construção matemática rigorosa: o **processo de Wiener** $W = (W_t)_{t \geq 0}$ é a idealização contínua do passeio aleatório. Hoje, é considerado o processo mais central e estudado da probabilidade ([verbete Wiener process](https://en.wikipedia.org/wiki/Wiener_process), cuja leitura recomendamos fortemente).

**Definição (padrão).** Um processo real $(W_t)_{t \geq 0}$ é um **movimento browniano padrão** quando:

1. $W_0 = 0$;
2. $W$ tem incrementos independentes: $W_t - W_s$ independente de $(W_u)_{u \leq s}$ para $s < t$;
3. Incrementos estacionários gaussianos: $W_t - W_s \sim N(0,\, t - s)$ — a variância cresce linearmente com a duração do intervalo;
4. As trajetórias $t \mapsto W_t(\omega)$ são contínuas (quase certamente).

Observações importantes:

- $\mathrm{Var}(W_t) = t$: a distribuição de $W_t$ "esparrama" conforme $t$ cresce; não há estado estacionário (é um processo não estacionário!) mesmo sendo suas leis de incremento estacionárias.
- Covariância: $\mathrm{Cov}(W_s, W_t) = \min(s, t)$ — a dependência entre dois instantes é proporcional ao menor deles.
- É **martingala** quadrado-integrável ($E[W_t \mid \mathcal{F}_s] = W_s$) e satisfaz a propriedade de Markov forte: a partir de qualquer tempo de parada, o processo reinicia como novo BM.
- Sua **variação quadrática** vale $[W]_t = t$ — este fato discreto-but-deep (a soma dos quadrados dos incrementos de uma partição tende ao limite $t$) é o motor do cálculo estocástico de Itô: sem ele, nada de lema de Itô, Black–Scholes etc. Comparações acessíveis em [discussão StackExchange sobre a variação quadrática do Wiener process](https://math.stackexchange.com/questions/3705345/computation-of-the-quadratic-variation-of-wiener-process) e [nota introdutória sobre variação quadrática em finanças](https://www.finance-tutoring.fr/the-relation-between-brownian-motion-and-quadratic-variation-in-layman%E2%80%99s-terms%E2%80%A6?mobile=1).

**Como simular.** Como os incrementos são Gaussianos com variância igual à duração, basta gerar incrementos normais independentes:

```{r}
# --- Caminho do movimento browniano padrão -------------------------------
set.seed(11)
T  <- 5                      # horizonte
dt <- 0.005                  # passo temporal
n  <- floor(T / dt)
inc <- rnorm(n, sd = sqrt(dt))       # incremento ~ N(0, dt)
W <- c(0, cumsum(inc))               # acumula (com W_0 = 0)
tgrid <- seq(0, T, length.out = n + 1)

plot(tgrid, W, type = "l", col = "firebrick", lwd = 1.6,
     xlab = "tempo", ylab = "W(t)", main = "Movimento browniano padrão (um caminho)")
abline(h = 0, col = "gray60")
lines(tgrid, 2 * sqrt(tgrid), col = "gray50", lty = 2)    # faixa ±2 desvio-padrão
lines(tgrid, -2 * sqrt(tgrid), col = "gray50", lty = 2)
```

Verifique as propriedades com seu próprio código: rode o bloco acima 2000 vezes com `replicate()` e confira (i) $E[W_t] \approx 0$; (ii) $\mathrm{Var}(W_t) \approx t$ para $t = 1, 2, 4$; (iii) a fração de caminhos que cruzam o nível $+2$ antes do tempo $t=4$ versus a previsão $P(W_4 > 2) = 1 - \Phi(1) \approx 0{,}159$ (usando reflexão no máximo). 

**Do passeio ao BM (escala).** Defina $W^{(n)}_t := \frac{1}{\sqrt{n}} S_{\lfloor nt \rfloor}$, onde $S$ é passeio aleatório simétrico. Então, conforme $n \to \infty$, os caminhos de $W^{(n)}$ convergem em lei para o movimento browniano (teorema de Donsker). Isso é mais que uma curiosidade: significa que **BM é o limite universal de processos de soma de passos pequenos**, explicando sua ubiquidade. Informalmente, é o TLC em dimensão funcional.

**Variações usadas na prática:**

- **BM com deriva e volatilidade:** $X_t = \mu t + \sigma W_t$ — tendência determinística + flutuação Browniana; aproxima preços de commodities, taxas de câmbio, etc.
- **Movimento geométrico browniano (GBM):** $S_t = S_0\exp\bigl((\mu - \tfrac{\sigma^2}{2})t + \sigma W_t\bigr)$ — solução da EDO estocástica $dS = \mu S\,dt + \sigma S\,dW$ (lema de Itô). O termo $-\sigma^2/2$ é a correção de Jensen/Itô: o valor esperado de $\exp(\sigma W_t)$ exige compensar a convexidade; sem ele, "modelaria errado". Este é o modelo por trás de Black–Scholes e de incontáveis simuladores de Monte Carlo em finanças.
- **BM fracionário** (parâmetro de Hurst $H$): generalização com memória longa; $H=1/2$ é o caso padrão. Usado em modelos de tráfego de rede, hidrologia e rough volatility. (Menção apenas; requer ferramentas específicas.)

```{r}
# --- GBM: um caminho de preço com drift mu e vol sigma -------------------
set.seed(21)
S0 <- 100; mu <- 0.08; sig <- 0.25; T <- 2; nstep <- 1000
dt <- T / nstep
inc <- rnorm(nstep, mean = (mu - sig^2 / 2) * dt, sd = sig * sqrt(dt))
S <- S0 * exp(cumsum(c(0, inc)))
plot(seq(0, T, length.out = nstep + 1), S, type = "l", col = "navy", lwd = 2,
     xlab = "tempo", ylab = "preço", main = "Caminho simulado de GBM")
abline(h = S0, lty = 2, col = "gray50")
abline(h = S0 * exp(mu * T), lty = 3, col = "firebrick")   # tendência teórica
```

**Aplicações típicas:** difusão molecular e nanopartículas; precificação e risco em finanças; ruído em sistemas dinâmicos; campos aleatórios Bayesianos; filtragem não-linear; física estatística.

### 4.5 Processos de renovação (trocar o relógio da exponencial)

**Motivação.** No Poisson, os tempos entre eventos são exponenciais — memória sem memória. Mas em muitos sistemas reais isso é irrealista: lâmpadas têm desgaste (falhas mais prováveis após certo uso), máquinas têm manutenções programadas. Os **processos de renovação** permitem distribuições arbitrárias (mas ainda i.i.d.) para os tempos entre eventos ([Renewal theory](https://en.wikipedia.org/wiki/Renewal_theory); notas práticas em [notas de renovação da Columbia](http://www.columbia.edu/~ks20/4106-18-Fall/Notes-Renewal-Theory.pdf)).

**Definição.** Dados i.i.d. $\tau_1, \tau_2, \dots$ com função de distribuição $F$ (concentrada em $(0,\infty)$), definimos os tempos de renovação $T_n = \sum_{k=1}^n \tau_k$ e o processo de contagem $N(t) = \max\{n : T_n \leq t\}$ — o número de renovações completadas até $t$. O caso $F = \mathrm{Exp}(\lambda)$ recupera o Poisson; escolhas determinísticas dão processos periódicos; gamas/Erlang dão "memória tipo-envelhecimento" (hipoexponencial).

Resultados centrais (declarados informalmente, todos demonstráveis nas fontes citadas):

- **Teorema elementar da renovação:** $E[N(t)] \sim t / E[\tau]$ quando $t \to \infty$ (fora casos aritméticos degenerados). Ou seja, assintoticamente vale a regra intuitiva "eventos por unidade de tempo = velocidade inversa da duração típica do ciclo".
- **Paradoxo da inspeção:** olhando o processo num instante fixo (digamos, agora), o intervalo em curso NÃO tem distribuição $F$ típica — é enviesado para intervalos longos. Formalmente, a distribuição do comprimento observado é proporcional a $z\,F(dz)/E[\tau]$. Consequência prática famosa: se você chega num ponto de ônibus aleatório e os intervalos entre ônibus são irregulares, seu tempo médio de espera excede metade do intervalo típico — chegar "no meio de um intervalo longo" é mais provável. Da mesma forma, a "vida residual" de uma lâmpada que já funcionou $a$ horas depende de $a$ de maneira explícita:
$$P(\text{resíduo} > x \mid \text{idade} > a) = \frac{1 - F(a + x)}{1 - F(a)} \quad \text{(contínuo)},$$
e a idade de um item em operação típicamente excede a vida média pela quantidade ligada à variabilidade de $F$.
- **Renovação com recompensas:** se cada ciclo rende "recompensa" $R_k$ (custo de manutenção, energia gerada etc.), o ganho médio de longo prazo é $\frac{E[R_1]}{E[\tau_1]}$ — princípio que sustenta otimização de políticas de substituição em engenharia de confiabilidade.

```{r}
# --- Renovação: Exp (Poisson) vs. Gama (não-exponencial) -----------------
set.seed(99)
horizonte <- 2000
Ns_em_uma_repeticao <- function(taus, t){ sum(cumsum(taus) <= t) }

# caso exponencial (média 1 -> taxa efetiva 1)
r_exp <- replicate(50, {
  taus <- rexp(10000, rate = 1)
  Ns_em_uma_repeticao(taus, horizonte)
})
# caso gama(2,1) (média 2 -> taxa efetiva 1/2)
r_gam <- replicate(50, {
  taus <- rgamma(10000, shape = 2, rate = 1)
  Ns_em_uma_repeticao(taus, horizonte)
})

c(EXP_media = mean(r_exp) / horizonte,
  EXP_teoria = 1,          # 1/E[tau] com tau ~ Exp(1)
  GAM_media = mean(r_gam) / horizonte,
  GAM_teoria = 0.5)        # 1/E[tau] com tau ~ Gama(2,1)
```

O exercício é interpretar a saída e comparar com o teorema elementar: a razão $E[N(t)]/t$ deve se aproximar de $1/E[\tau]$. Note também que o erro de Monte Carlo aqui deve ser reportado (seção 7) — use intervalos de confiança via `t.test` sobre as replicações.

**Aplicações:** análise de disponibilidade de sistemas reparáveis; política de manutenção preventiva; seguros (frequência de sinistros não-Poissonianos); economia (substituição de equipamentos); redes (distribuição de interrupções).

### 4.6 Séries temporais (observar uma única trajetória)

Até aqui, em geral tínhamos "muitas réplicas do mesmo experimento". Em séries temporais temos o cenário oposto: **uma única trajetória de comprimento $n$** (a chuva diária de São Paulo por 30 anos, a demanda mensal por 200 meses) — e queremos aprender o processo subjacente a partir dela. Séries temporais são o ramo da modelagem estocástica dedicado a dados indexados pelo tempo.

**Estacionariedade (em sentido amplo).** Um processo $\{Z_t\}_{t \in T}$ é **fracamente estacionário** quando $E[Z_t] = \mu$ (constante), $E[Z_t^2] < \infty$ e a covariância $\gamma(s,t) = \mathrm{Cov}(Z_s, Z_t)$ depende apenas do atraso $|s - t|$ ([Notas de processos estocásticos IME-USP](https://www.ime.usp.br/~lane/home/MAE5870/STprocestv2.pdf)). Se a distribuição conjunta completa é invariante por translação de tempo, dizemos **estritamente estacionária** (mais restritiva). **Ergodicidade**, grosso modo, é a propriedade que permite "ler" estatísticas do processo inteiro a partir de uma trajetória longa — o que legitima estimativas de médias/correlações a partir de dados históricos ([slides UFSC](https://www.inf.ufsc.br/~marcelo.menezes.reis/Processos01.pdf)).

**Função de autocorrelação (FAC).** Para uma série estacionária, define-se

$$\rho(k) = \frac{\gamma(k)}{\gamma(0)}, \qquad k = 0, 1, 2, \dots,$$

onde $\gamma(k) = \mathrm{Cov}(Z_t, Z_{t+k})$. Ela mede quanto a série "se lembra" de si própria após $k$ períodos: tipicamente, decai a zero (devagar, se memória longa; rápido, se memória curta), ou oscila com picos sazonais. É o instrumento central de diagnóstico em séries temporais — e, na prática, estimada a partir de uma amostra finita (`acf()` no R).

**Modelos ARMA(p,q).** Um processo **ARMA**$(p,q)$ combina autorregressão (dependência direta dos $p$ valores passados) e médias móveis (dependência dos $q$ choques passados):

$$Z_t = \phi_1 Z_{t-1} + \cdots + \phi_p Z_{t-p} + \theta_1 a_{t-1} + \cdots + \theta_q a_{t-q} + a_t,$$

com $\{a_t\}$ i.i.d. (geralmente Normal, ruído branco). Hipóteses clássicas de **causalidade/invertibilidade** garantem representação em termos do passado e do futuro e estão ligadas às raízes das equações características ficarem fora do círculo unitário. Casos notáveis: AR(1) com $\phi_1 \neq 0$ gera correlação que decai geometricamente; MA(1) tem FAC que corta abruptamente após lag 1; ARMA(1,1) mistura comportamentos. Para séries com estrutura sazonal, usa-se SARIMA$(p,d,q)(P,D,Q)_m$, onde $m$ é a periodicidade (12 para mensal, 7 para semanal).

**Metodologia de Box-Jenkins (o ciclo prático):**

1. **Identificar** ordem/estrutura a partir de FAC/FACP e transformações necessárias (diferenciar para remover tendência/estacionarizar; log para estabilizar variância);
2. **Estimar** coeficientes por mínimos quadrados/máxima verossimilhança (pacote `forecast::auto.arima` automatiza seleção heurística; sempre conferir resíduos!);
3. **Verificar** — teste Ljung-Box na FAC dos resíduos, gráficos de normalidade e heterocedasticidade; se insuficiente, voltar ao passo 1;
4. **Prever** — equações de previsão com intervalos de confiança crescentes com o horizonte (a incerteza aumenta com o lead time; nunca extrapole longe demais com confiança exagerada).

**Quando a estacionariedade falha.** Muitos fenômenos econômicos/climáticos não são estacionários: a variância cresce, o nível muda, a estrutura se altera. Dois perigos clássicos:

- **Regressão espúria:** regredir duas séries não estacionárias pode produzir R² alto e coeficientes significantes mesmo sem qualquer relação causal real ([MathWorks — Time Series Regression IV: Spurious Regression](https://www.mathworks.com/help/econ/time-series-regression-iv-spurious-regression.html); [estudo sobre não-estacionaridade](https://mpra.ub.uni-muenchen.de/27926/1/Stationarity_of_time_series_and_the_problem_of_spurious_)).
- **Raiz unitária / integração:** testes (Dickey-Fuller aumentado etc.) verificam se a série possui componente estocástico permanente que exige diferenciação ([notas CGSpace/ CGIAR sobre pitfalls de raiz unitária](https://cgspace.cgiar.org/server/api/core/bitstreams/7aa92438-0578-4751-b149-6a8f8f12dfb8/content)).

A receita prática: sempre examine FAC, tente transformar/diferenciar, valide out-of-sample — e documente decisões, porque pré-tratamento escolhido post-hoc introduz viés de pós-escolha (ligação com a seção 7).

```{r}
# --- Branco, cinza e vermelho: três exemplos sintéticos de espectro -------
set.seed(1)
n <- 256
branco <- rnorm(n)                                  # ruído branco: FAC ~ delta
cinza  <- arima.sim(list(order = c(1, 0, 0), ar = 0.7), n = n)   # AR(1): decaimento
serie  <- cumsum(rnorm(n))                          # passeio aleatório: não-estacionária

par(mfrow = c(2, 3))
plot(branco, main = "ruído branco"); plot(cinza, main = "AR(1) phi=0.7"); plot(serie, main = "passeio aleatório (cumsum)")
acf(branco, lag.max = 30, main = "FAC branco")
acf(cinza,  lag.max = 30, main = "FAC AR(1)")
acf(serie,  lag.max = 30, main = "FAC passeio")
```

Compare os painéis: o decaimento da FAC identifica a memória do AR(1); a FAC do ruído branco não mostra padrão sistemático; a FAC do passeio quase não decai — sinal clássico de não-estacionaridade.

### Exercícios da seção 4 (continuação)

5. Num processo de Poisson de taxa $\lambda = 2$/min, qual é $P(N(5) = 0)$? E o tempo esperado até a primeira chegada? Simule e compare.
6. Mostre empiricamente (simulação) que $\mathrm{Var}(W_t)/t \to 1$ para o movimento browniano. Qual é o efeito de discretizar com $\Delta t$ grande?
7. Explique o paradoxo da inspeção com o exemplo de um trem que passa a cada 10 minutos (determinístico) versus a cada $\mathrm{Exp}(1/10)$ min (Poisson). Se você chega num instante aleatório, qual é o tempo médio de espera nos dois casos?
8. Para uma série de preços de ação (dados reais baixados de alguma fonte), faça a FAC e discuta: estacionária em nível? em retornos? que transformação seria defensável?

*Soluções comentadas no Apêndice D.*


# 5. Simulação estocástica em R: do básico ao reproduzível

Simular processos estocásticos é uma das melhores formas de *entendê-los* — e uma habilidade central para quem trabalha com dados, risco ou modelagem científica. Nesta seção reunimos o ferramental prático: geração de números aleatórios, ferramentas base do R, os ecossistemas especializados (séries temporais, espacial, Bayesiano) e o que separa um script profissional de um resultado irreprodutível. O script `demo_processos_estocasticos.R` acompanha este documento e pode ser rodado ponta a ponta.

## 5.1 Aleatoriedade controlada: `set.seed()`

O R gera números "pseudoaleatórios" por algoritmos determinísticos — o gerador padrão é o **Mersenne Twister MT19937**. "Pseudo" significa: se você sabe a semente, reproduz exatamente a mesma sequência. Por isso, **toda simulação séria começa com `set.seed()`**:

```{r}
set.seed(2026)             # fixa a semente da sequência inteira
amostra <- sample(1:6, size = 100, replace = TRUE)   # 100 lançamentos de dado
media_sim <- mean(amostra)

# Repita a linha acima quantas vezes quiser: mesmo resultado sempre.
```

Boas práticas e sutilezas importantes:

- Fixe **uma vez** no topo do script; resetar dentro de loops quando não há motivo claro apenas esconde variabilidade.
- A reprodutibilidade exige também **mesma versão do R e dos pacotes** usados ([Random Seeds and Reproducibility](https://medium.com/data-science/random-seeds-and-reproducibility-933da79446e3); [Choosing the "Correct" Seed for Reproducible Research, StackExchange](https://stats.stackexchange.com/questions/335936/choosing-the-correct-seed-for-reproducible-research-results)). Use `sessionInfo()` ou, melhor ainda, um snapshot de bibliotecas (`renv::snapshot()`).
- Sementes são arbitrárias — nenhuma é "melhor" a priori; escolha uma e mantenha-a honestamente registrada ([Stack Overflow: reasons for using set.seed](https://stackoverflow.com/questions/13605271/reasons-for-using-the-set-seed-function)).
- Em **paralelismo**, cada trabalhador precisa de sua própria semente/stream — nunca replique a mesma semente em todos os núcleos, ou suas estimativas serão correlacionadas e você perderá precisão sem perceber ([Setting a seed in R when using parallel simulation](https://www.r-bloggers.com/2018/07/%F0%9F%8C%B1-setting-a-seed-in-r-when-using-parallel-simulation/)). Uma solução simples: derive sementes por índice, `seed_base + i`, e registre-as; soluções mais robustas usam fluxos independentes (método de L'Ecuyer), disponíveis em pacotes como `future`.

## 5.2 Ferramentas básicas de simulação

A maior parte deste guia usa apenas funções fundamentais do R, mostradas abaixo. Dominar esses cinco minutos de sintaxe libera qualquer modelo.

**(a) `sample()` e `replicate()`: Monte Carlo na veia.**

```{r}
# Média empírica da soma de 100 dados — converge para 350 pela LGN
sim <- replicate(10_000, sum(sample(1:6, size = 100, replace = TRUE)))
hist(sim, breaks = 30, freq = FALSE, main = "Soma de 100 dados: 10.000 simulações")
abline(v = 350, col = "red", lwd = 2)
c(media = mean(sim), desvio_padrao = sd(sim))
```

Cuidado comum: `sample(x)` sem reposição é **permutação**; `sample(x, replace = TRUE)` é **bootstrap**. Confundir os dois invalida muitos resultados clássicos (teste permutacional × bootstrap são coisas diferentes!).

**(b) Famílias de geradores:** `runif`, `rnorm`, `rexp`, `rbeta`, `rgamma`, `rpois`, `rt` etc. — todas seguem o padrão `r<nome>` (gerar) / `d<nome>` (densidade) / `p<nome>` (acumulada) / `q<nome>` (quantil). Aprenda o padrão; vale para tudo.

**(c) `arima.sim()`: séries temporais sintéticas.** Gera séries de modelos ARIMA/SARIMA especificados — excelente para testar metodologia antes de aplicar em dados reais:

```{r}
set.seed(7)
serie <- arima.sim(list(order = c(2, 1, 1), ar = c(0.8, -0.2), ma = 0.3), n = 500)
plot(serie, main = "Série sintética ARIMA(2,1,1) — 500 pontos")
```

**(d) `simulate()`: distribuição preditiva de modelos ajustados.** Para um modelo estatístico já ajustado (`lm`, `glm`, ...), gera novas respostas plausíveis sob os pressupostos estimados (por padrão, reamostragem dos resíduos):

```{r}
# Exemplo com modelo linear sobre grade de cenários climáticos:
dados_ficticios <- data.frame(
  temp = rep(c(18, 22, 26), each = 10),
  rain = rep(c(900, 1300), each = 15),
  yield = c(rnorm(15, mean = 5, sd = 1), rnorm(15, mean = 7, sd = 1))
)
modelo <- lm(yield ~ temp + rain, data = dados_ficticios)
cenarios <- expand.grid(temp = seq(15, 30, by = 5), rain = c(800, 1200, 1600))
pred_sim <- simulate(modelo, nsim = 200, newdata = cenarios)
summary(pred_sim)     # intervalos preditivos por cenário
```

**(e) Vetorização e velocidade.** Prefira operações vetoriais/matríciais a loops explícitos; use `apply`/`lapply`/`sapply`/`purrr::map`; para cargas pesadas considere `Rcpp` (C++ compilado) ou paralelismo. Regra prática: primeiro faça funcionar corretamente; depois otimize (mas sempre reporte custo computacional e tempo de execução no seu protocolo).

## 5.3 Séries temporais: classes e infraestrutura

Para análise de séries temporais no R, a trilha moderna é mapeada pela [CRAN Task View: Time Series Analysis](https://cran.r-project.org/view=TimeSeries) (mantida por Rob Hyndman e Rebecca Killick, atualizada continuamente). Mapa resumido de pacotes e funções essenciais:

| Tarefa | Pacotes principais | Funções-chave |
|---|---|---|
| Classes e manipulação | `ts`, `xts`, `anytime`, `tsibble`, `fabletools` | `ts()`, `as_tsibble()`, `feasts`, `modeltime` |
| Decomposição/suavização | `stats` (base), `forecast` | `decompose()`, `stl()`, `ses()`, `holt()`, `ets()` |
| Modelos ARIMA/SARIMA | `forecast`, `tseries`, `sarima` | `auto.arima()`, `arima()`, `sarima.fit()` |
| Diagnóstico | `forecast`, `tseries`, `tsdiag` | `acf()`, `pacf()`, `Box.test()`, `adf.test()`, `tsdiag()` |
| Previsão com incerteza | `forecast`, `fable` | `forecast()`, `predict()`, `interval = TRUE` |
| Dados reais (exemplos) | `datasets`, `fpp3`, `dadosBR` | `load`, séries econômicas/britânicas |

Dica canônica: trabalhe num único fluxo consistente — carregue → transforme → ajuste → valide out-of-sample → documente. Evite misturar três APIs diferentes no mesmo projeto; o esforço migra para a lógica. Para séries **espacialmente referenciadas** ("clima em estações"), adicione os componentes espaciais da próxima subseção.

## 5.4 Geoespacial: padrões pontuais, campos e krigagem

Duas famílias de problemas aparecem constantemente em climatologia, epidemiologia, ecologia e geociências:

**(a) Padrões pontuais** (eventos localizados: epicentros, árvores doentes, crimes, raios): use **`spatstat`** — representações `ppp`, medidas de intensidade, distâncias ao vizinho mais próximo, função $K$ de Ripley e envelopes de completude espacial aleatória ([Getting started with spatstat](https://cran.r-project.org/web/packages/spatstat/vignettes/getstart.pdf); capítulo dedicado em [Spatial Statistics for Data Science](https://www.paulamoraga.com/book-spatial/the-spatstat-package.html)). Fluxo típico:

```{r}
# if (!requireNamespace("spatstat", quietly = TRUE)) install.packages("spatstat")
library(spatstat)
win <- owin(xrange = c(0, 1), yrange = c(0, 2))       # janela de observação
X <- runifpoint(n = 200, win = win)                   # 200 pontos uniformes
plot(X, main = "Padrão pontual simulado (uniforme)")
intensity(X)                 # intensidade média estimada
nn <- nndist(X)              # distâncias ao vizinho mais próximo
hist(nn$dist, main = "Distância ao vizinho mais próximo")
```

Testes formais (`Kest()`, `envelope()`) decidem se os dados exibem agregação, regularidade ou aleatoriedade completa — fundamentos para modelagem de risco espacial.

**(b) Campos contínuos** (concentração de poluentes, altitude, precipitação): aqui entra a **geoestatística** — o variograma mede como a semelhança entre locais decai com a distância; a krigagem usa essa estrutura para interpolar/prever em locais não amostrados, com intervalos de incerteza. A implementação clássica é o pacote **`gstat`** ([tutorial gstat/INPE](http://www.dpi.inpe.br/gilberto/tutorials/software/R-contrib/gstat.pdf)); a teoria acessível está em [capítulo Kriging do livro de Moraga](https://www.paulamoraga.com/book-spatial/kriging.html), que inclui o exemplo clássico dos dados de zinco na planície aluvial do rio Meuse:

```{r}
# if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat")
library(gstat)
data(meuse)                       # 155 amostras de metais pesados (Meuse, NL)
var_raw <- variogram(log(zn) ~ 1, locations = ~x + y, data = meuse)
plot(var_raw, main = "Variograma experimental (log zinco)")

var_fit <- fit.variogram(var_raw, vgm(psill = 0.5, model = "Exp", range = 500))
print(var_fit)                    # efeito pepita, patamar, alcance

grid <- meuse.grid                # grade de predição disponível com o pacote sp
zk <- krige(log(zn) ~ 1, locations = ~x + y, data = meuse, newdata = grid, model = var_fit)
spplot(zk["var1.pred"], main = "Krigagem ordinária: log zinco previsto")
```

Interprete: o **alcance** indica até onde pontos estão correlacionados; o **efeito pepita** captura microestrutura/erro de medição; a krigagem atribui pesos maiores a vizinhos próximos e informa a variância de predição (mapa de incerteza) — nada menos que um modelo probabilístico completo do campo.

O ecossistema vai longe: `sf` (dados vetoriais), `terra`/`raster` (campos raster), `spdep`/`INLA` (modelos Bayesianos espaciais), além da tarefa temporal-espaço combinada ([CRAN SpatioTemporal task view](https://cran.r-project.org/view=SpatioTemporal)).

## 5.5 Amostragem MCMC e inferência Bayesiana com `rstan`

Quando o modelo tem muitos parâmetros e dados pequenos, ou quando queremos incerteza completa (distribuições posteriores inteiras), entram os métodos **MCMC** (*Markov Chain Monte Carlo*). O Stan é hoje o motor mais usado, e o R dialoga com ele via **`rstan`** ([vinheta oficial RStan](https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html)):

``` stan
// modelo.stan — Normal-Normal hierárquico bem simples
data {
  int<lower=1> N;            // nº de observações
  array[N] real y;           // dados
}
parameters {
  real mu;                   // média populacional
  real<lower=0> sigma;       // desvio-padrão
}
model {
  sigma ~ exponential(1);    // prior suave
  mu ~ normal(0, 10);
  y ~ normal(mu, sigma);     // verossimilhança
}
```

```{r}
# if (!requireNamespace("rstan", quietly = TRUE)) install.packages("rstan")
# library(rstan)
# fit <- stan(file = "modelo.stan", data = list(N = length(y), y = y),
#             iter = 2000, warmup = 1000, chains = 4, seed = 2026)
# print(fit)                      # resumo posterior: média, dp, ICs 95%
# summary(fit)$summary[, "Rhat"]  # diagnóstico de convergência (deve ser ~1)
```

Pontos de atenção (detalhados na própria vinheta e na [documentação de MCMC do Stan](https://mc-stan.org/docs/reference-manual/mcmc.html)): o sampler padrão (**NUTS**) faz adaptação de passo durante o *warmup*; cadeias múltiplas devem convergir entre si; monitore `Rhat` e tamanho efetivo de amostra (**ess**), não apenas o número bruto de iterações; divergências indicam geometria problemática do posterior (reespecifique o modelo ou reparametrize). Alternativas amigáveis: `coda` para MCMC clássico, `brms` para modelos de regressão Bayesianos "de fórmulas", `INLA` para inferência Bayesiana aproximada escalável em geoestatística.

## 5.6 Redução de variância (para fazer mais com menos)

Erro Monte Carlo cai tipicamente na ordem $1/\sqrt{N}$ — dobrar a precisão requer quadruplicar as simulações ([Monte Carlo Convergence and Standard Error, Risk Hub](https://riskhub.org/blogs/monte-carlo-convergence-and-standard-error)). Três técnicas clássicas reduzem a constante multiplicativa (sem sacrificar validade):

- **Variáveis antitéticas:** simule pares opostos (por exemplo, $Z$ e $-Z$) — caminhos negativamente correlacionados cancelam parte da flutuação. Em precificação de opções, funciona notavelmente bem porque payoffs monotônicos ficam altamente correlacionados.
- **Variáveis de controle:** use uma quantidade correlacionada com o alvo, cujo valor esperado é conhecido, para corrigir o estimador (um análogo de regressão para Monte Carlo).
- **Amostragem estratificada:** divida o espaço amostral em estratos e aloque simulações proporcionalmente; garante cobertura de caudas e reduz variância ([Variance Reduction Techniques for Efficient Monte Carlo](https://mbrenndoerfer.com/writing/variance-reduction-monte-carlo-simulation)).

## 5.7 Protocolo reprodutível (checklist)

Um checklist que usamos na escrita deste documento e recomendamos em todo projeto estocástico:

1. **Uma semente, registrada** — `set.seed()` no início; anotar valor no texto/resultados.
2. **Script único e ordenado** — carregar → simular → analisar → plotar; nada interativo dependente de estado global.
3. **Dados imutáveis** — guardar dados brutos com hash/versão; derivar cópias por código; nunca editar planilha "à mão".
4. **Versões registradas** — `sessionInfo()`, `renv::snapshot()`, Docker/containers quando possível; documentar versões do sistema.
5. **Erros reportados** — intervalos de Monte Carlo, número de replicações, tempo de CPU; nunca reporte um único caminho como "o resultado".
6. **Sanity checks** — verifique casos-limite analíticos conhecidos (probabilidade 0/1, esperanças teóricas) antes de confiar no pipeline.
7. **Rodada final limpa** — apague artefatos antigos, rode o script do zero, confira saídas iguais (bits idênticos para os mesmos seeds e versões).

Isso conversa diretamente com boas práticas gerais de ciência aberta e reprodutibilidade ([Reprodutibilidade na ciência — USP](https://cienciaaberta.usp.br/reprodutibilidade-na-ciencia/); [Guias práticos para uma pesquisa mais reprodutível](https://www.reprodutibilidade.org/post/guias-praticos-para-uma-pesquisa-mais-reprodutivel)).

### Exercícios da seção 5

1. Estime $\pi$ usando `replicate()` de médias de `sample(1:6, n, replace=TRUE)` para `n = 10, 100, 1000` e mostre graficamente a concentração em torno de 3,5. Comente a relação com a LGN.
2. Reproduza a série `arima.sim(list(order=c(1,0,0), ar=0.9), n=1000)` e estime a FAC. Qual é o comportamento teórico esperado? Teste com outro valor de `ar` negativo.
3. Faça um estudo pequeno de redução de variância: compare o erro-padrão de `mean(exp(Z))` com `Z ~ N(0,1)` estimado com 10.000 simulações comuns vs. 5.000 pares antitéticos. Qual é o ganho observado?
4. Escreva um plano de reprodução para os experimentos deste capítulo (sementes, versões, gráficos). O que você faria diferente se fosse publicar?

*Soluções comentadas no Apêndice D.*


# 6. Aplicações guiadas (casos completos em R)

Esta seção reúne quatro aplicações completas — finanças, epidemias, filas e clima — no mesmo espírito dos exemplos anteriores: problema → modelo → código → interpretação. Todo o código está consolidado no script `demo_processos_estocasticos.R`. **Aviso:** os parâmetros são didáticos; nada aqui é recomendação de investimento.

## 6.A Finanças: precificação de opções por Monte Carlo

**Problema.** Uma opção de compra europeia dá o direito de comprar um ativo por $K$ (preço de exercício) ao vencimento $T$. Se soubermos a dinâmica do preço (volatilidade $\sigma$, taxa livre de risco $r$), qual é o "preço justo" hoje? Em 1973, Black & Scholes deram uma fórmula fechada para o caso lognormal — mas Monte Carlo ensina a mecânica geral que funciona quando a fórmula não existe (opções americanas, exóticas, múltiplos ativos).

**Modelo.** Sob a medida risk-neutral (mundo onde o retorno esperado do ativo é $r$), o preço segue movimento geométrico browniano:

$$dS_t = r S_t\,dt + \sigma S_t\,dW_t \quad\Longrightarrow\quad S_T = S_0\exp\Bigl((r - \tfrac{\sigma^2}{2})T + \sigma W_T\Bigr),$$

e o preço da opção é o valor presente do payoff descontado: $C_0 = e^{-rT}\,E\bigl[\max(S_T - K, 0)\bigr]$. A expectativa é estimada por média sobre caminhos simulados — simples assim.

```{r}
# --- Precificação Monte Carlo de opção de compra europeia ----------------
set.seed(123)

precificar_call_mc <- function(S0 = 100, K = 100, r = 0.05, sigma = 0.20,
                               T = 1, npath = 20000, nstep = 500){
  dt <- T / nstep
  caminho <- matrix(S0, nrow = npath, ncol = nstep + 1)
  for (k in 1:nstep){
    inc <- rnorm(npath, mean = (r - sigma^2 / 2) * dt, sd = sigma * sqrt(dt))
    caminho[, k + 1] <- caminho[, k] * exp(inc)
  }
  payoff <- pmax(caminho[, nstep + 1] - K, 0)
  preco  <- exp(-r * T) * mean(payoff)
  erro   <- exp(-r * T) * sd(payoff) / sqrt(npath)
  c(preco = preco, erro_padrao = erro, npath = npath)
}

# Fórmula de Black-Scholes para conferência (call europeia)
bs_call <- function(S, K, r, sigma, T){
  d1 <- (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
}

mc  <- precificar_call_mc()
th  <- bs_call(100, 100, 0.05, 0.20, 1)
c(mc, black_scholes = th)
# Esperado: MC ≈ 10,4–10,5 (erro-padrão ~0,2) e BS ≈ 10,44
```

**Extensão opcional: variáveis antitéticas.** Onde o payoff é suave e monotônico, a versão antitética converte o ruído com muito mais eficiência: simule metade dos incrementos, derive a outra metade trocando o sinal, e calcule a média dos dois payoffs aos pares:

```{r}
set.seed(124)
call_antiteticas <- function(S0 = 100, K = 100, r = 0.05, sigma = 0.20,
                             T = 1, npares = 10000, nstep = 500){
  dt <- T / nstep
  payoff_pos <- numeric(npares); payoff_neg <- numeric(npares)
  for (k in 1:nstep){
    inc <- rnorm(npares, mean = (r - sigma^2 / 2) * dt, sd = sigma * sqrt(dt))
    pos <- cumsum(inc); neg <- -pos
    payoff_pos <- pmax(S0 * exp(pos[nstep]) - K, 0)
    payoff_neg <- pmax(S0 * exp(neg[nstep]) - K, 0)
  }
  par <- (payoff_pos + payoff_neg) / 2      # cada par forma uma média
  c(preco = exp(-r * T) * mean(par), erro_padrao = exp(-r * T) * sd(par) / sqrt(npares))
}
call_antiteticas()
```

Discuta: compare tempo × precisão dos três métodos (MC bruto, MC antitético, BS). Por que a diferença antitética diminui com passos pequenos? (O payoff de cada caminho fica altamente correlacionado.) Quando BS não existir, o MC ainda roda; essa é a razão de seu sucesso industrial ([Black–Scholes vs Monte Carlo na prática](https://www.acturtle.com/blog/call-option-valuation-black-scholes-vs-monte-carlo); [exemplo implementado em R, StackExchange Quant](https://quant.stackexchange.com/questions/50168/monte-carlo-option-pricing-with-r)).

**Sanity checks obrigatórios:** (i) sem aleatoriedade ($\sigma = 0$), o resultado deve ser $\max(S_0 e^{rT} - K, 0)$ descontado; (ii) $C_0 \le S_0$ sempre; (iii) comparando com BS, a diferença deve cair como $1/\sqrt{n_{\text{paths}}}$.

## 6.B Epidemias: o modelo SIR estocástico

**Problema.** Uma doença contagiosa se espalha numa população de $N$ indivíduos. Como prever se haverá surto, quão grande será, e quanto tempo dura? O modelo **SIR** divide a população em três compartimentos — Suscetíveis (S), Infectados (I) e Recuperados/imunes (R) — com fluxos entre eles.

**Equações determinísticas (fluxo médio):**

$$\frac{dS}{dt} = -\beta S I, \qquad \frac{dI}{dt} = \beta S I - \gamma I, \qquad \frac{dR}{dt} = \gamma I,$$

com taxas $\beta$ (transmissão) e $\gamma$ (recuperação), e fração inicial de infectados pequena. Um número adimensional decide o destino do sistema: o **número básico de reprodução**

$$R_0 = \frac{\beta N}{\gamma},$$

que conta quantos casos secundários, em média, um infectado gera numa população totalmente suscetível. Se $R_0 < 1$: o surto morre quase certamente; se $R_0 > 1$: cresce rápido, tipicamente até consumir grande parte dos suscetíveis antes de desacelerar (a famosa curva achatada). Nas versões estocásticas — cadeias de Markov de tempo contínuo (CTMC) ou equações diferenciais estocásticas —, tudo isso vale com probabilidade: há chance de extinção mesmo com $R_0 > 1$ (mais relevante em populações pequenas), e o limiar fica borrado perto de $R_0 = 1$ ([compartmental models](https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)); [tutorial do King sobre simulação estocástica, com código R no pacote pomp](https://kingaa.github.io/short-course/stochsim/stochsim.html); aplicação brasileira ao sarampo em [Lau et al., TEMA/SBMAC](https://tema.sbmac.org.br/tema/article/download/1670/1173)).

```{r}
# --- SIR determinístico x estocástico (Gillespie simplificado) -----------
SIR_det <- function(N = 1000, beta = 0.002, gamma = 0.0015,
                    I0 = 5, Tmax = 200, dt = 0.1){
  n <- floor(Tmax / dt)
  S <- N - I0; I <- I0; R <- 0
  traj <- matrix(NA, n + 1, 3)
  traj[1, ] <- c(S, I, R)
  for (k in 1:n){
    dI <- (beta * S * I / N - gamma * I) * dt
    dS <- -beta * S * I / N * dt
    dR <- gamma * I * dt
    S <- max(S + dS, 0); I <- max(I + dI, 0); R <- R + dR
    traj[k + 1, ] <- c(S, I, R)
  }
  data.frame(t = seq(0, Tmax, length.out = n + 1), S = traj[,1], I = traj[,2], R = traj[,3])
}

SIR_est <- function(N = 1000, beta = 0.002, gamma = 0.0015, I0 = 5, Tmax = 200){
  # CTMC: eventos de infecção e recuperação com taxas instantâneas
  S <- N - I0; I <- I0; R <- 0
  saida <- list(); t <- 0
  while (t < Tmax && I > 0){
    tempos <- rexp(2, rate = c(beta * S * I / N, gamma * I))   # próximos eventos
    tipo  <- which.min(tempos); dt_ <- min(tempos)
    if ((t + dt_) > Tmax) break
    t <- t + dt_
    if (tipo == 1) { S <- S - 1; I <- I + 1 } else { I <- I - 1; R <- R + 1 }
    saida[[length(saida) + 1]] <- c(t, S, I, R)
  }
  as.data.frame(do.call(rbind, c(list(c(0, N - I0, I0, 0)), saida)))
}

det <- SIR_det(); st <- SIR_est()
plot(det$t, det$I, type = "l", col = "darkred", lwd = 2, ylab = "infectados",
     main = "SIR: fluxo determinístico (linha cheia) vs. trajetória estocástica (pontilhada)")
lines(st$t, st$I, col = "steelblue", lty = 2)
legend("topright", legend = c("I determinístico", "I estocástico"),
       col = c("darkred", "steelblue"), lty = c(1, 2))
```

**Pontos para discussão:** (i) rode com $R_0 = 1{,}33$ (`beta=0.002, gamma=0.0015` em $N=1000$) e observe o crescimento; (ii) reduza para $R_0 < 1$ e note a extinção quase certa; (iii) em populações pequenas ($N=20$), mesmo com $R_0 > 1$, há extinguições por acaso; (iv) o pico estocástico oscila em torno do determinístico — a diferença é o "ruído demográfico"; (v) modelos realistas acrescentam exposição sazonal, idade, rede de contatos etc.; o ferramental bayesiano permite calibrar $\beta$ e $\gamma$ a partir de dados (via `pomp`/Stan — seção 5.5).

## 6.C Filas: o M/M/1 (por que você espera tanto?)

**Problema.** Clientes chegam a uma central de atendimento segundo um processo de Poisson de taxa $\lambda$ (por minuto) e são atendidos com tempos exponenciais de média $1/\mu$. Qual é o tempo médio de espera? Quantos clientes ficam na fila? Quando o sistema explode?

**Teoria (fórmulas fechadas!).** Para o modelo M/M/1 com estabilidade $\rho := \lambda/\mu < 1$ ([M/M/1 queue — Wikipedia](https://en.wikipedia.org/wiki/M/M/1_queue)):

$$P_0 = 1 - \rho \quad\text{(prob. de sistema vazio)}, \qquad L_q = \frac{\rho^2}{1-\rho} \quad\text{(clientes médios na fila)}$$

$$W_q = \frac{\rho^2}{\mu(1-\rho)} \quad\text{(tempo médio de espera na fila)}, \qquad W = W_q + \frac{1}{\mu}, \qquad L = \lambda W.$$

Note a estrutura: o tempo no sistema = espera + serviço; a lei de Little ($L = \lambda W$) aparece como identidade consistência (aproveite para verificá-la empiricamente!). Note também: quando $\rho \to 1$ (utilização alta), os tempos de espera crescem explosivamente — daí a sensação cotidiana de "quase ninguém aguenta mais esperar". O estado $\rho \ge 1$ faz a fila crescer indefinidamente: instável.

```{r}
# --- M/M/1: teoria versus simulação --------------------------------------
mm1_teorico <- function(lambda, mu){
  rho <- lambda / mu
  if (rho >= 1) stop("rho >= 1: fila instável!")
  Wq <- rho^2 / (mu * (1 - rho))
  W  <- Wq + 1 / mu
  c(rho = rho, P_vazio = 1 - rho, Lq = lambda * Wq, Wq = Wq, W = W, L = lambda * W)
}

mm1_sim <- function(lambda, mu, nclientes = 5000, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  chegadas <- cumsum(rexp(nclientes, rate = lambda))     # instantes de chegada
  servicos <- rexp(nclientes, rate = mu)                 # tempo de atendimento
  inicio_atend <- numeric(nclientes); fim_atend <- numeric(nclientes)
  for (i in 1:nclientes){
    inicio_atend[i] <- max(chegadas[i], if (i == 1) 0 else fim_atend[i - 1])
    fim_atend[i]    <- inicio_atend[i] + servicos[i]
  }
  espera <- (inicio_atend - chegadas)[fim_atend > 0]
  c(Wq_media = mean(espera), W_media = mean(fim_atend - chegadas))
}

mm1_teorico(lambda = 3, mu = 5)     # espera na fila 0,18 min; Lq = 0,54
mm1_sim(lambda = 3, mu = 5, seed = 1)
```

Explore depois:

- O efeito de subir $\rho$ para 0,9: os tempos multiplicam por ~4,5;
- M/M/c (vários servidores): queda dramática de espera ([guia aplicado M/M/1 vs M/M/c](https://metricgate.com/blogs/queueing-theory-mm1-mmc-models/));
- Distribuição não-exponencial do serviço (M/G/1 — fórmula de Pollaczek–Khinchine exige outros momentos; use a notação de Kendall e discuta aproximações).

## 6.D Clima: chuvas de São Paulo — estudo de caso com dados reais

**Dados.** Utilizamos o banco climático do projeto (reanálise ERA5 via Open-Meteo; tabela mensal por município, cobrindo todos os 5.571 municípios brasileiros). Consultamos precipitação mensal acumulada para São Paulo (capital; id IBGE 3550308), resultando em séries anuais 1995–2025 e mensais recentes. Os arquivos foram salvos em `dados/chuva_anual_sp.csv`, `dados/chuva_mensal_sp_recente.csv` e `dados/chuva_climatologia_mensal_sp_1991_2020.csv`. **Limitação importante:** reanálise climática combina modelo e observações interpoladas, com resolução ~9–25 km, o que suaviza extremos locais (chuva de tempestade intensa pode ser subestimada vs. estação física do INMET; temperaturas urbanas podem aparecer ~1 °C menores). Use-a para análises estruturais e ensino; para extremos críticos, valide com medições de superfície.

**Panorama descritivo (números calculados na consulta):**

- Chuva anual 1995–2025: média **1254,6 mm**, mediana 1266,4 mm, desvio-padrão **186,4 mm**; mínimo **864,4 mm (2014)**, máximo **1646,9 mm (1995)**. A série sugere flutuações decadais — médias por década: 1332 mm (1995–2004), 1250 mm (2005–2014) e 1203 mm (2015–2024); últimos 10 anos completos: 1158 mm. Autocorrelação lag-1 estimada: 0,36 (persistência moderada ano a ano).
- Climatologia mensal 1991–2020: verão chuvoso (jan 222,3 mm; fev 165,4; mar 139,0; dez 169,7) e inverno seco (ago 35,8 mm; jun 44,5; jul 48,9), típico do regime de monção subtropical. Máximas históricas em janeiro chegaram a 374,5 mm e mínimas a 0,7 mm em julho — a amplitude é enorme.
- Recente (mensal): dezembro/2024 teve **293,8 mm**; já junho/2024 registrou apenas 2,8 mm. Em 2025 choveu 1038 mm no total; 2026 acumula **1216,1 mm até outubro**, com setembro/outubro muito chuvosos (222 e 221 mm), sinal de recuperação após meses secos.

> **Figura de referência (dados reais):** ![](doc/figs/chuva_sp_anual_e_climatologia.png)

```{r}
# --- Estudo de caso: chuva anual de São Paulo (dados reais) --------------
anos  <- 1995:2025
chuva <- read.csv("dados/chuva_anual_sp.csv")$chuva_mm

x <- ts(chuva, start = 1995, frequency = 1)
plot(x, main = "Chuva anual acumulada — São Paulo (1995–2025)", ylab = "mm")
abline(h = mean(chuva), col = "red", lty = 2)              # média histórica

par(mfrow = c(1, 2))
acf(x, lag.max = 10, main = "FAC da chuva anual (memória moderada)")
pacf(x, lag.max = 10, main = "FACP")

set.seed(2026)
fit <- arima(x, order = c(1, 1, 0))                        # ARIMA(1,1,0) parcimonioso
print(fit)                                                 # phi e constante
checkresiduals(fit)                                        # diagnóstico (Figura)
```

```{r}
# --- Climatologia mensal e previsão one-step -----------------------------
clim <- read.csv("dados/chuva_climatologia_mensal_sp_1991_2020.csv")
barplot(clim$chuva_media_hist, names.arg = month.name[1:12],
        col = "cornflowerblue", ylab = "mm", main = "Climatologia mensal (1991–2020)")

recente <- read.csv("dados/chuva_mensal_sp_recente.csv")
y <- ts(recente$chuva_mm, start = c(2024, 1), frequency = 12)
plot(y, main = "Chuva mensal recente (2024–2026)")
abline(h = 120, col = "gray50", lty = 2)                   # referência visual
```

Interpretar com honestidade: o padrão sazonal é robusto; as variações decadais existem mas exigem testes formais (Mann-Kendall, Sen's slope, análise de changepoints) antes de afirmar "tendência"; amostras de ~30 pontos têm pouca potência para separar sinal de ruído vermelho. Isso nos leva diretamente às armadilhas da próxima seção.

### Exercícios da seção 6

1. No exemplo de opções, teste $\sigma = 0{,}4$ e veja como o erro Monte Carlo cresce; explique por que. Depois compare com o preço BS correspondente.
2. No SIR, estime a probabilidade de o surto se extinguir começando com $I_0 = 1$ para diferentes $R_0$ (faça 500 simulações para cada valor) — trace o gráfico da probabilidade de extinção versus $R_0$.
3. Para o M/M/1 com $\lambda = 4$ e $\mu = 6$, confira analiticamente que $L = \lambda W$ e meça o tempo de simulação necessário para atingir erro-padrão pequeno.
4. No estudo climático, faça a FAC das primeiras diferenças $(X_t - X_{t-1})$ e verifique se algum padrão sazonal remanescente aparece; proponha um SARIMA$(0,1,0)(0,1,1)_{12}$ para os dados mensais e discuta seus resíduos.

*Soluções comentadas no Apêndice D.*


# 7. Armadilhas, erros comuns e controvérsias

Nenhuma técnica é neutra. Modelos estocásticos são ferramentas de comunicação sobre incerteza — e como toda ferramenta, podem ser mal-usados. Esta seção documenta armadilhas conhecidas, com fontes para leitura posterior. **Leia antes de tirar conclusões de qualquer análise estocástica.**

## 7.1 Sementes, gerações de números aleatórios e falsa reprodutibilidade

Já dissemos na seção 5 que fixar semente é obrigatório; aqui vão as sutilezas que causam retrabalho:

- **Re-seedar dentro de loop sem intenção** cria dependências artificiais e elimina variação legítima. Cada `set.seed()` define um ponto de partida da sequência; se você re-executa o mesmo bloco com a mesma semente, obtém o mesmo resultado (bom para debug) — mas se fizer isso dentro de uma rotina de estimação, está congelando os dados falsos.
- **Paralelismo ingênuo:** dar a mesma semente a todos os trabalhadores produz replicações correlacionadas — o erro aparente cai, mas o erro real não ([Setting a seed in R when using parallel simulation](https://www.r-bloggers.com/2018/07/%F0%9F%8C%B1-setting-a-seed-in-r-when-using-parallel-simulation/)). Use fluxos derivados corretamente (`future::future.seed`, `L'Ecuyer`), ou derive sementes por índice e registre todas.
- **Versões diferentes:** trocar versão do R ou de pacote pode mudar o comportamento de funções e até o algoritmo de amostragem subjacente. Reprodutibilidade verdadeira exige versionamento de ambiente ([Random Seeds and Reproducibility](https://medium.com/data-science/random-seeds-and-reproducibility-933da79446e3); [Choosing the "Correct" Seed — Cross Validated](https://stats.stackexchange.com/questions/335936/choosing-the-correct-seed-for-reproducible-research-results)).
- **"Escolher a semente pelo resultado"** é uma forma sutil de p-hacking: testar várias sementes e reportar a que dá o resultado desejado infla falsos positivos. A semente deve ser escolhida *antes* — idealmente sorteada ao ar livre, literalmente (dado de seis lados para gerar dígitos, como sugerido nas discussões sobre escolha "honesta" de sementes).
- **Pseudo-aleatoriedade ≠ aleatoriedade física:** geradores são determinísticos; para uso criptográfico, não servem (não é nosso caso, mas vale saber). Os resultados "exatos" de um script nada têm de "prova" do fenômeno modelado — só das consequências do modelo sob aquele conjunto específico de hipóteses.

## 7.2 P-valores, "significância" e o jardim dos caminhos

O p-valor mede **quão surpreendente** é o resultado observado sob uma hipótese específica (tipicamente "nada acontecendo"). Ele **não** mede: (i) a probabilidade de a hipótese nula ser verdadeira; (ii) o tamanho do efeito; (iii) a utilidade prática do achado. Confundir essas coisas está na origem de décadas de crise estatística ([The Statistical Crisis in Science — American Scientist](https://www.americanscientist.org/article/the-statistical-crisis-in-science)).

Uma armadilha particularmente insidiosa é o **"jardim dos caminhos que se bifurcam"** (garden of forking paths) descrito por Andrew Gelman e Eric Loken ([texto completo](https://sites.stat.columbia.edu/gelman/research/unpublished/p_hacking.pdf)): mesmo sem nenhuma "pesca" consciente de resultados, cada decisão de análise — qual variável incluir, transformar ou não, excluir outliers, cortar a amostra em subgrupos — cria ramificações; quando as decisões são contingentes aos dados, dezenas de análises defensáveis podem coexistir, e basta uma delas resultar $p < 0{,}05$ para aparecer conclusão. O problema não requer má-fé: apenas flexibilidade razoável somada à ausência de pré-especificação.

Complementos práticos: simulações sistemáticas documentam que "p-hacking" amplo (testar muitas especificações e reportar o melhor) infla enormemente taxas de falso positivo ([Big little lies: a compendium and simulation of p-hacking — Royal Society Open Science](https://royalsocietypublishing.org/rsos/article/10/2/220346/92017/Big-little-lies-a-compendium-and-simulation-of-p)), e propostas de reforma incluem baixar o limiar convencional de 0,05 para 0,005 ([Impact of redefining statistical significance on p-hacking — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11098386/)). Você pode simular o fenômeno você mesmo:

```{r}
# --- Demonstração: em dados puramente aleatórios, algum teste fica "significativo"
set.seed(99)
nrep <- 1000; nobs <- 120

p_min <- replicate(nrep, {
  x <- rnorm(nobs)                       # nenhum sinal estruturado
  ps <- c(
    t.test(x)$p.value,                     # média global
    t.test(x[seq_len(nobs) %% 2 == 1])$p.value,  # subgrupo "ímpar"
    t.test(x[seq_len(nobs) %% 2 == 0])$p.value,  # subgrupo "par"
    cor.test(x, lag(x))$p.value            # autocorrelação lag-1
  )
  min(ps)                               # melhor (menor) p entre 4 tentativas
})

hist(p_min, breaks = 20, col = "gray80",
     main = "Menor p-valor entre 4 testes — dados 100% aleatórios")
abline(v = 0.05, col = "red", lwd = 2)
mean(p_min < 0.05)   # tipicamente > 0,5: quase sempre "alguém" passa
```

A moral: **um único resultado significativo não é evidência contra o acaso se o procedimento de análise tinha liberdade escondida**. Antídotos conhecidos: pré-registro do plano de análise; análise de sensibilidade ("multiverse"); reportar tamanhos de efeito com intervalos de confiança; distinguir claramente exploração de confirmação; e nunca usar termos como "prova" para $p < 0{,}05$.

Erro correlato: **barras de erro enganadoras**. Cumming e colaboradores ([Error bars in experimental biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC2064100/); [Robust misinterpretation of confidence intervals — PMC](https://link.springer.com/article/10.3758/s13423-013-0572-3)) documentaram que barras de erro e intervalos de confiança são rotineiramente mal interpretados — inclusive a regra errada de "se as barras de dois grupos se sobrepõem, não há diferença significativa". Sempre explicite o tipo de barra (dp, erro-padrão, IC 95%, intervalo preditivo), o nível, e a unidade; e prefira apresentar intervalos estimados a depender apenas de regras de ouro opacas.

## 7.3 Estacionariedade: assumir sem verificar

Grande parte da teoria clássica (ARIMA, krigagem ordinária, médias móveis) pressupõe estacionariedade de alguma forma. Dados reais frequentemente violam: tendências, sazonalidades instáveis, mudanças estruturais (mudanças de regime no clima, crises econômicas), ciclos com períodos variáveis. As consequências são graves e bem documentadas:

- **Regressão espúria:** duas séries ambas com tendência produzem correlações altas "significativas" mesmo sem relação causal ([MathWorks — Spurious Regression](https://www.mathworks.com/help/econ/time-series-regression-iv-spurious-regression.html); [Baumöhl — Stationarity of time series and the problem of spurious regression](https://mpra.ub.uni-muenchen.de/27926/1/Stationarity_of_time_series_and_the_problem_of_spurious_)).
- **Raiz unitária e inferência distorcida:** testes usuais (t, F) deixam de ter distribuição padrão quando componentes permanentes existem; o caminho correto envolve diferenciação ou modelagem de longo prazo ([Unit Root Tests: Common Pitfalls and Best Practices — CGIAR](https://cgspace.cgiar.org/server/api/core/bitstreams/7aa92438-0578-4751-b149-6a8f8f12dfb8/content)).
- **Confundir os tipos de não-estacionaridade:** série "com tendência" e série "integrada" exigem tratamentos diferentes; remover tendência de algo integrado não resolve (e vice-versa). A distinção exige testes formais + teoria, não apenas gráficos ([Time Series Nonstationarity — bookdown R for Political Scientists](https://bookdown.org/ccolonescu/RPoE4/time-series-nonstationarity.html)).
- **Transformações milagrosas:** diferenciar demais ("over-differencing") induz memória longa artificial e complica a interpretação; log ajuda em variâncias crescentes, mas muda a escala de interpretação; filtrar/desenhar séries "limpas" esconde escolhas subjetivas. Documente cada transformação com justificativa.

## 7.4 Erros típicos em simulações Monte Carlo

Além das armadilhas gerais de simulação (mitos e erros clássicos catalogados há décadas — [INFORMS/WSC 2001 paper](https://informs-sim.org/wsc01papers/006.PDF); [On the Assessment of Monte Carlo Error — Am. Statistician](https://pmc.ncbi.nlm.nih.gov/articles/PMC3337209/)):

1. **Reportar uma única trajetória** como "o resultado". Simulação estocástica é coleção de trajetórias; um caminho individual é anedótico.
2. **Ignorar o erro Monte Carlo**: todo número estimado vem com incerteza $\propto 1/\sqrt{N}$; reporte SEM/ICs e verifique convergência ([Monte Carlo Convergence and Standard Error — Risk Hub](https://riskhub.org/blogs/monte-carlo-convergence-and-standard-error)). O famoso padrão de "quadruplicar amostras para dobrar precisão" decorre disso.
3. **Subestimar variância por amostrar errado:** reamostragens dependentes, blocos insuficientes em MCMC, ou variáveis antitéticas usadas sem contabilizar a correlação levam a intervalos enganosamente estreitos.
4. **Calibrar e simular com a mesma semente/família de caminhos** pode enviesar comparações entre políticas ("common random numbers" ajudam a comparar alternativas pareadas, mas devem ser usados com desenho adequado).
5. **Distribuições de entrada erradas:** caudas leves onde há saltos (falhas, crashes, epidemics superdisseminadores) subestimam riscos extremos ([Monte Carlo Simulation — practical guide](https://mcpanalytics.ai/articles/monte-carlo-simulation-practical-guide-for-data-driven-decisions)). Faça análise de sensabilidade às distribuições e use testes de adequação (KS, QQ-plots).
6. **Achar que mais iterações resolvem tudo:** se o modelo está mal especificado, velocidade computacional apenas acelera a resposta errada. Validação face à realidade vem antes de otimização numérica.

## 7.5 Confusões conceituais frequentes

**(a) Determinístico × estocástico × caótico.** Um sistema pode ser determinístico e ainda parecer aleatório (caos determinístico: pequenas diferenças de condição inicial amplificam exponencialmente; exemplo clássico: mapa logístico $x_{n+1} = 4x_n(1-x_n)$). Nem toda irregularidade é ruído; nem toda aleatoriedade observada vem de randomness fundamental. Em finanças, o passeio aleatório foi proposto como modelo de preços ([Fama — Random Walks in Stock-Market Prices](https://www.chicagobooth.edu/~/media/34F68FFD9CC04EF1A76901F6C61C0A76.PDF)), mas a própria hipótese de mercado eficiente segue debatida quanto a previsibilidade e bolhas ([Investopedia — Random Walk Theory](https://www.investopedia.com/terms/r/randomwalktheory.asp); [Samuelson vs Fama on the Efficient Market Hypothesis — Œconomia](https://journals.openedition.org/oeconomia/5300)). Lição metodológica: distinguir **modelo** de **mecanismo**; o passeio aleatório é um excelente modelo *nulo* (baseline), não a verdade última do mercado.

**(b) Poisson × interchegadas exponenciais × renovação.** Usar Poisson quando as inter-chegadas têm memória (envelhecimento/devido a desgaste) viola a propriedade sem memória e subestima risco; usar renovação simples para eventos mutuamente excitantes (uma falha provoca outra) também quebra. Diagnóstico: ajuste distribuicional das inter-chegadas (QQ-plots), testes de explosão/rejeição da constância de taxa.

**(c) Independência onde há dependência.** Muitas simplificações (binomial vs. hipergeométrico, M/M/1 vs. G/G/k, amostragem com reposição) assumem independência condicional correta. Testar suposições de independência (run tests, autocorrelação, mutual information) antes de aplicar modelos i.i.d.

**(d) confundir tempo discreto e contínuo.** Equações em diferenças vs. EDOs/EDEs; discretizações grosseiras mudam dinâmica (por ex., no SIR com passo temporal grande, populações negativas aparecem — cuidado!). Sempre reportar o esquema de discretização adotado.

**(e) confundir "processo" e "trajetória".** O processo é a lei probabilística; a trajetória é uma realização. Estatísticas de uma única trajetória podem não representar a lei (problema de ergodicidade quebrada em tempos longos — sistemas com dois regimes que não se misturam, vidros girando lentamente etc.). Antes de calcular médias temporais longas, verifique ergodicidade.

## 7.6 Controvérsias ativas (para onde olhar quando avançar)

1. **Mercados eficientes vs. previsibilidade:** debate vivo desde Fama/Samuelson; hoje inclui fatores de risco sistemáticos, *momentum*, e crítica pós-2008 sobre modelos de risco gaussiano. Consenso parcial: mercados são "quase" eficientes em janelas curtas, com desvios persistentes em horizontes longos.
2. **Reforma da inferência estatística:** substituição de limiares dicotômicos ($p < 0{,}05$) por estimações com incerteza; Bayesian hierarchical models; pré-registro; ciência "open science framework". Críticas de ambos os lados continuam em publicações.
3. **Volatilidade e roughness em finanças:** modelos de vol fracionária e "rough volatility" contestam a suavidade dos modelos clássicos; implicações para precificação e gestão de risco ainda em disputa.
4. **Longa memória e Hurst em hidroclimatologia:** séries de vazão/chuva podem exibir comportamento fractal — discussão sobre métodos de detecção de tendência e periodicidades naturais.
5. **Inferência causal em séries temporais:** diferenças-em-diferenças, controles sintéticos e descontinuidades são objetos ativos de pesquisa; escolha de defasagens e testes de robustez seguem debatidos.
6. **IA e quantificação de incerteza:** redes neurais com estimativas Bayesianas (*Bayesian neural networks*), conformal prediction, e validação de modelos generativos são fronteiras com muito trabalho pela frente.

## 7.7 Checklist final de qualidade

Um mini-exame rápido antes de aceitar qualquer análise estocástica:

- [ ] Modelo e hipóteses explicitados; testes de adequação feitos (resíduos, FAC, distribuições)?
- [ ] Transformações de dados justificadas e documentadas?
- [ ] Estacionariedade verificada (ou modelada) com testes + teoria?
- [ ] Dependências (temporal/espaçal) consideradas?
- [ ] Incerteza propagada: erro Monte Carlo reportado, ICs preditivos, análise de sensibilidade?
- [ ] Validação out-of-sample / backtesting realizada?
- [ ] Reprodução completa: seeds, versões, pipeline imutável?
- [ ] Linguagem: "resultado consistente com X sob as hipóteses Y", nunca "prova"?

Se cinco ou mais caixas estiverem verdes, ótimo. Se menos, trate os resultados como exploratórios — com humildade pública.

### Exercícios da seção 7

1. Replique a demonstração do jardim dos caminhos com oito testes (adicione testes em quartis, médias movidas, etc.) e mostre que $\Pr(\min p < 0{,}05)$ cresce. Qual seria a correção de Bonferroni nesse caso e por que ela é conservadora?
2. Gere duas séries integradas independentes (passeios aleatórios) de comprimento 200 e faça a regressão $Y_t = \alpha + \beta X_t + \varepsilon_t$; repita para 100 pares e faça um histograma dos R². Compare com o resultado esperado de séries relacionadas (mesmo caminho comum). Discuta "regressão espúria".
3. Num modelo AR(1) com $\phi = 0{,}95$, estime por OLS a autocorrelação lag-1 em amostras de tamanho 50; discuta viés em amostras pequenas versus a teoria assintótica.
4. Explique por que as barras de erro do tipo "erro-padrão da média" ficam mais estreitas que os intervalos de predição, e quando cada uma é apropriada.

*Soluções comentadas no Apêndice D.*


# 8. Para continuar estudando

Selecionamos materiais a partir das frentes investigadas nesta documentação: cursos brasileiros com ementas e notas disponíveis publicamente, tutoriais internacionais de referência, recursos para aprender **R** aplicado e caminhos de estudo estruturados. Todos os itens abaixo foram localizados e verificados na web (links diretos).

## 8.1 Trilha em português

### Cursos universitários com material aberto

| Recurso | Conteúdo | Link |
|---|---|---|
| **IME-Unicamp — ME501 Processos Estocásticos** (graduação) | Ementa completa: cadeias de Markov discretas, classificação de estados, medida estacionária, processos de Poisson, nascimento-morte, noções de filas e simulação; bibliografia clássica (Feller, Hoel-Port-Stone, Karlin-Taylor) | [Página oficial](https://www.ime.unicamp.br/graduacao/disciplina/me501-processos-estocasticos) |
| **Unicamp — MI414 Introdução aos Processos Estocásticos** (pós) | Elementos de processos estocásticos, definições, classificações, confiabilidade e estabilidade | [Página oficial](https://www.ime.unicamp.br/pos-graduacao/disciplina/mi414-introducao-aos-processos-estocasticos) |
| **USP — MAE0312** | Cadeia de Markov, comportamento assintótico, processo de Poisson | [Registro JupiterWeb](https://uspdigital.usp.br/jupiterweb/obterDisciplina?sgldis=mae0312) |
| **USP — ME5870 / Slides "Processos Estocásticos"** (IME-USP) | Transparências cobrindo estacionariedade fraca, funções de correlação, processos Gaussianos | [Slides STprocestv2.pdf](https://www.ime.usp.br/~lane/home/MAE5870/STprocestv2.pdf) |
| **Notas de aula do Prof. Antonio Galves (IME-USP)** | Aplicações de processos estocásticos; estimação de cadeias de Markov via máxima verossimilhança, entropia | [PDF](https://www.ime.usp.br/~galves/cursos/NotasdeAulaAplProcEst.pdf) |
| **UFMG — Conceitos Básicos de Processos Estocásticos** (RTE-01/2011, Atuncar) | Notas didáticas sobre definição, exemplos típicos e construção de processos; referência amigável ao iniciante | [PDF](https://www.est.ufmg.br/portal/wp-content/uploads/2023/01/RTE-01-2011.pdf) |
| **UFSC — Apostila de Processos Estocásticos** | Visão geral com ênfase em classificação (inclui estacionariedade e ergodicidade) | [PDF](https://www.inf.ufsc.br/~marcelo.menezes.reis/Processos01.pdf) |
| **UFOP — Capítulo 6: Processos Estocásticos e Cadeias de Markov** | Material introdutório para engenharia de produção, com exercícios | [PDF](https://www.professor.ufop.br/sites/default/files/magno/files/capitulo_6_-_processos_estocasticos_e_cadeias_de_markov_1.pdf) |
| **UFC — TI0112 Processos Estocásticos** | Plano de ensino: definições, processos estacionários e ergódicos, Poisson, cadeias de Markov, filas | [PDF](https://engcomp.ufc.br/wp-content/uploads/2020/08/ti0112-processos-estocasticos.pdf) |
| **FGV — Processos Estocásticos (EMAP)** | Programa com bibliografia obrigatória (Ross) | [Página](https://emap.fgv.br/disciplina/processos-estocasticos-1) |
| **PUC-Rio — capítulo “3 Processos Estocásticos”** | Definição e exemplos iniciais em nível aplicado | [PDF](https://www.maxwell.vrac.puc-rio.br/15321/15321_4.PDF) |
| **CEPEL/SEE — Definições acerca de Processos Estocásticos** | Explicação didática contrastando processos determinísticos e estocásticos | [Página](https://see.cepel.br/manual/anatem/recursos/serie/estocastico.html) |
| **SBMAC — Monografias (ex.: livro 81)** | Coleção de monografias em matemática aplicada, alguns volumes sobre processos estocásticos | [PDF](https://www.sbmac.org.br/wp-content/uploads/2022/08/livro_81.pdf) |

Para quem quiser ir além, o livro **"Introdução aos Processos Estocásticos" de Daniel Miranda (UFABC)** — justamente o que usamos como fio condutor deste documento — cobre o assunto com uma abordagem moderna, diagramas interativos e foco em intuição ([capítulo 1](https://danielmiranda.prof.ufabc.edu.br/livros/processos/capitulo-01-processos-estocasticos.html)). A obra faz parte do movimento de livros-texto abertos e interativos no Brasil.

### Aprender R para dados e séries temporais

| Recurso | Conteúdo | Link |
|---|---|---|
| **Curso-R — Introdução à Programação com R** | Curso inicial brasileiro, slides e material prático | [Slides](https://curso-r.github.io/main-intro-programacao/slides/01-introducao-ao-curso.pdf) |
| **Escola Virtual de Governo — Análise de Dados em Linguagem R** | Do básico ao machine learning, com certificado | [Página do curso](https://www.escolavirtual.gov.br/curso/325) |
| **Guiader (GitHub)** | Guia colaborativo de trilhas, cursos, livros e canais em português | [Repositório](https://github.com/arthurspk/guiader) |
| **Alura — Guia Completo da Linguagem R** | Panorama de pacotes, séries temporais e econometria em R | [Artigo](https://www.alura.com.br/artigos/linguagem-r) |
| **Análise de Séries Temporais usando R (analisemacro)** | Curso aplicado de séries temporais econométricas | [Página](https://analisemacro.com.br/curso/analise-de-series-temporais-usando-r/) |
| **Introdução às Séries Temporais com R (Medium)** | Artigo prático com código e datasets | [Texto](https://medium.com/@flaviagaia/introdu%C3%A7%C3%A3o-%C3%A0s-s%C3%A9ries-temporais-com-a-linguagem-r-37157435cfbd) |
| **Séries Temporais em R (rpubs)** | Tutorial comentado | [Publicação](https://rpubs.com/RobsonOL/series_temporais_r) |
| **Epidemiology Handbook em PT — capítulo de séries temporais** | Detecção de surtos e séries temporais aplicadas à saúde | [Capítulo](https://www.epirhandbook.com/pt/new_pages/time_series.pt.html) |
| **Livro aberto "Spatial Statistics for Data Science" (Moraga et al.)** | Capítulos sobre krigagem com gstat e análise de padrões pontuais com spatstat | [Krigagem](https://www.paulamoraga.com/book-spatial/kriging.html) · [spatstat](https://www.paulamoraga.com/book-spatial/the-spatstat-package.html) |

### Clássicos em inglês recomendados pelas ementas brasileiras

- **William Feller**, *An Introduction to Probability Theory and Its Applications*, vol. 1 e 2 (Wiley) — referência histórica; cita-se frequentemente que a maioria dos processos aqui tratados já aparece nele.
- **Karlin & Taylor**, *A First Course in Stochastic Processes* (Academic Press) — texto usado nos cursos da Unicamp; cobre cadeias, renovais e aplicações.
- **Sheldon Ross**, *Introduction to Probability Models* — adotado pela FGV; muitos exemplos acessíveis com R.
- **Hoel, Port & Stone**, *Introduction to Stochastic Processes* (Houghton-Mifflin) — citado nas ementas do IMECC.

## 8.2 Tutoriais internacionais de referência (seleção com leitura feita)

- **MIT OpenCourseWare — 6.262 Discrete Stochastic Processes**: capítulos sobre Poisson e renovação, com desenvolvimento axiomático e problemas resolvidos ([cap. 2 — Poisson](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/3a19ce0e02d0008877351bfa24f3716a_MIT6_262S11_chap02.pdf); [cap. 4 — Renovação](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/931ffa0940899c27f34b71ad64fd2bb0_MIT6_262S11_chap04.pdf)). Excelente segundo passo depois das notas brasileiras.
- **Leiden University — Spieksma, *An Introduction to Stochastic Processes in Continuous Time***: notas europeias concisas sobre processos de Lévy, martingalas e cálculos estocásticos ([PDF](https://www.math.leidenuniv.nl/~spieksma/colleges/sp-master/LN.pdf)).
- **CMU Statistics 754 — notas de Shalizi**: construções via distribuições de dimensões finitas, tempos de parada e filtros ([lecture 2](https://www.stat.cmu.edu/~cshalizi/754/notes/lecture-02.pdf); [lecture 6](https://www.stat.cmu.edu/~cshalizi/754/notes/lecture-06.pdf)).
- **randomservices.org — Filtrações e Tempos de Parada**: exposição rigorosa com demonstrações completas; fonte desta seção 3 ([capítulo](https://www.randomservices.org/random/prob/Stop.html)); também [introdução à renovação](https://www.randomservices.org/random/renewal/Introduction.html).
- **UChicago (Lalley) — notas de Brownian motion e martingalas**: duas fontes usadas nas seções 4 e 3 ([BM](https://galton.uchicago.edu/~lalley/Courses/313/BrownianMotionCurrent.pdf); [martingalas](http://galton.uchicago.edu/~lalley/Courses/385/Martingales.pdf)).
- **Purdue (Stindel) — Construction of the Wiener process**: roteiro construtivo passo a passo ([PDF](https://www.math.purdue.edu/~stindel/teaching/ma539/brownian-motion1.pdf)).
- **Dartmouth Math 100 — projetos finais** usados como apoio nas seções 4.1 e 3.2: [recorrência de Pólya (Mare)](https://math.dartmouth.edu/~pw/math100w13/mare.pdf); [teorema da parada opcional (LaLonde)](https://math.dartmouth.edu/~pw/math100w13/lalonde.pdf).
- **Jeremy Kun — Martingales and the Optional Stopping Theorem**: porta de entrada suave e divertida ([link](https://www.jeremykun.com/2014/03/03/martingales-and-the-optional-stopping-theorem/)).
- **Stanford/Stat — Perlin et al., GARCH tutorial com R** ([SciELO](https://www.scielo.br/j/rac/a/FyNM3rThZN7wQgMKshwLHjk/?lang=en)) para finanças quantitativas em R.
- **King et al. — Simulation-based Inference** ([material do short course](https://kingaa.github.io/short-course/stochsim/stochsim.html)): ponte entre modelagem mecanística e inferência Bayesiana.
- **CRAN Task Views** ([lista completa](https://cran.r-project.org/web/views/); [Time Series](https://cran.r-project.org/view=TimeSeries); [Espaço-Temporal](https://cran.r-project.org/view=SpatioTemporal); [Finance](https://cran.r-project.org/view=Finance)) — mapas vivos do ecossistema R.

## 8.3 Comunidades e prática

- **Perguntas e respostas:** Stack Overflow (tags `r`, `time-series`, `stochastic-process`, `markov-chain`) e Cross Validated (stats.stackexchange.com) — busque antes de perguntar; quase sempre existe resposta detalhada.
- **Listas oficiais:** R-help e R-sig (estatística), além de comunidades locais (grupos de usuários R no Telegram/Discord e capítulos acadêmicos brasileiros).
- **Sociedades científicas:** SBMAC (Matemática Aplicada e Computacional) e ABE (Associação Brasileira de Estatística) promovem eventos com trilhas de processos estocásticos e séries temporais.
- **Prática deliberada:** refaça todos os scripts deste documento mudando parâmetros; depois implemente um modelo novo *do zero* (por exemplo, o modelo SIR com sazonalidade, ou urnas de Pólya) e valide contra resultados conhecidos. Só então avance para artigos com resultados novos.

## 8.4 Roteiro sugerido de 12 semanas

| Semanas | Tópico | Entregável |
|---|---|---|
| 1–2 | Probabilidade essencial: RV, distribuições, esperança/variância, LLN, TCL | Simulação de jogos de dados; histogramas convergindo para normal |
| 3–4 | Passeio aleatório, recorrência, ruína do jogador | Gráficos de trajetórias 1D/2D/3D; relatório da probabilidade de ruína vs. teoria |
| 5–6 | Cadeias de Markov (discretas e contínuas), distribuições estacionárias | Matriz de transição de um problema real (clima/tráfego); validação πP = π |
| 7 | Processo de Poisson e renovação | Simulação de centro de chamadas; comparação com teoria |
| 8–9 | Movimento browniano, martingalas, cálculo estocástico discreto | Simulação GBM; precificação simples de opções |
| 10 | Séries temporais: FAC/PACF, ARIMA, testes de raiz unitária | Modelagem da chuva de São Paulo com backtest honesto |
| 11 | Simulação avançada: redução de variância, MCMC, paralelismo | Estudo Monte Carlo do método escolhido; relatório de desempenho |
| 12 | Projeto final + armadilhas | Relatório com checklist da seção 7.7 aplicado a si mesmo |

**Sugestão de projeto final:** prever o consumo elétrico semanal de uma cidade a partir de dados históricos, comparando três modelos (ARIMA/SARIMA, suavização exponencial, e um modelo mais complexo de sua escolha), com backtesting honesto, intervalos preditivos e discussão das limitações.

---

# Apêndice A — Glossário rápido

- **Processo estocástico:** família de variáveis aleatórias indexadas por algum conjunto (tempo, posição etc.).
- **Trajetória / caminho amostral:** função $t \mapsto X_t(\omega)$ obtida fixando um resultado $\omega$; um gráfico de simulação mostra uma trajetória.
- **Filtração $\{F_t\}$:** coleção crescente de σ-álgebras representando informação acumulada até cada instante.
- **Adaptado:** o valor $X_t$ é mensurável em relação à informação disponível até $t$.
- **Tempo de parada:** regra de decisão ("parar quando...") cuja ocorrência é decidível pelo histórico até cada momento; fundamental para apostas, testes sequenciais e finanças.
- **Martingala:** processo cuja expectativa condicional futura, dado o presente, é igual ao valor atual — modelo de jogo justo; pode ser sub- ou supermartingala conforme a expectativa sobe ou desce.
- **Teorema da parada opcional:** sob condições adequadas, a expectativa de uma martingala em um tempo de parada é igual à inicial (com requisitos de integrabilidade e cautela com tempos ilimitados).
- **Recorrente / transitivo:** um estado é recorrente se a caminhada volta a ele com probabilidade 1 (ou infinitas vezes); transitivo caso contrário. Em $d \geq 3$ dimensões, o passeio simétrico é transitivo (Pólya).
- **Distribuição estacionária $\pi$:** vetor de probabilidades tal que $\pi P = \pi$; distribuição limite de longo prazo para cadeias regulares.
- **Absorvente:** estado que, uma vez atingido, não se deixa mais; cadeia absorvente modela falhas, falências, extinções.
- **Taxa de transição (matriz $Q$, CTMC):** intensidades instantâneas de saltos entre estados; equações direta e inversa de Kolmogorov regem a evolução.
- **Processo de Poisson:** contagem de eventos raros e independentes com taxa constante $\lambda$; inter-chegadas i.i.d. exponenciais; propriedades: estacionariedade e independência dos incrementos.
- **Inter-chegadas:** tempos entre eventos consecutivos; no Poisson, $T_k - T_{k-1} \sim \mathrm{Exp}(\lambda)$ i.i.d.
- **Processo de renovação:** generaliza o Poisson permitindo qualquer distribuição das inter-chegadas (iid); modela desgaste, reparos, revisitas.
- **Paradoxo da inspeção:** o intervalo de renovação que contém um instante escolhido tende a ser maior que o intervalo típico ($E[\tau^2]/E[\tau] \geq E[\tau]$).
- **Vida residual / idade:** tempo remanescente de um ciclo em andamento; cresce com a idade no caso envelhecente.
- **Movimento browniano (Wiener):** processo Gaussiano de incrementos independentes, estacionários e médios zero; variância proporcional ao tempo; paths contínuos mas nunca diferenciáveis; protótipo universal (teorema de Donsker).
- **Variação quadrática:** soma dos quadrados dos incrementos; para BM, vale $[W]_t = t$ — pedra angular do cálculo estocástico e das opções financeiras.
- **Medida neutra ao risco / martingala descontada:** mundo em que preços descontados são martingalas; base da precificação arbitragem-free.
- **Equação de Black-Scholes / fórmula BS:** EDP/EDA que resulta em preço fechado para opções europeias lognormais; caso limite de outros métodos.
- **Modelo SIR:** compartimentos suscetível-infectado-recuperado; limiar de epidemia determinado por $R_0 = \beta N/\gamma$.
- **Número básico de reprodução $R_0$:** casos secundários esperados de um infectado isolado numa população inteiramente suscetível; limiar crítico 1.
- **CTMC:** cadeia de Markov a tempo contínuo com saltos instantâneos e taxas especificadas; base para modelos epidêmicos estocásticos e confiabilidade.
- **Simulação de Gillespie:** algoritmo exato (event-driven) para CTMC/reações químicas; alternativa rápida a Euler em sistemas rígidos.
- **SDE (equação diferencial estocástica):** equação com termo de ruído; integrações numéricas requerem cuidados especiais (convergência forte/fraca).
- **Lei dos Grandes Números:** médias amostrais convergem para a esperança verdadeira.
- **Teorema Central do Limite:** somas padronizadas de variáveis i.i.d. (ou dependentes com mixing) tendem a Normal; base para inferência e intervalos de confiança.
- **Autocorrelação (ACF/FAC):** correlação entre valores separados por defasagem $\tau$; mede memória temporal.
- **ARMA(p,q):** combinação autorregressiva e média móvel; modela séries estacionárias com memória curta.
- **SARIMA:** extensão sazonal do ARIMA para periodicidades conhecidas.
- **Raiz unitária / integração:** componente estocástico permanente; séries com tendência não-estacionária e validade de inferência modificada.
- **Diferenciação:** transformação $(1-B)^d X_t$ que remove memória longa e torna a série aproximadamente estacionária.
- **Erro-padrão (SE):** desvio-padrão do estimador; cai com $\sqrt{n}$; difere de erro-padrão Monte Carlo, que cai com $\sqrt{N_{\text{sims}}}$.
- **Monte Carlo:** estimativa de quantidades por reamostragem aleatória repetida; convergência $O(1/\sqrt{N})$.
- **Redução de variância:** técnicas (antitéticos, controles, estratificação) que diminuem o ruído de estimadores sem enviesá-los.
- **MCMC:** família de algoritmos (Metropolis, Gibbs, HMC/NUTS) que produzem amostras de distribuições-alvo via cadeias Markovianas ergódicas.
- **p-valor:** probabilidade de obter resultado extremo(s) assumindo hipótese nula; mede surpresa, não tamanho nem probabilidade da hipótese.
- **Jardim dos caminhos que se bifurcam:** multiplicidade de análises defensáveis contingentes aos dados que infla falsos positivos mesmo sem p-hacking explícito.
- **Regressão espúria:** associação aparente significativa entre séries não relacionadas, tipicamente por tendências comuns.
- **Ruído vermelho:** processo onde potências de baixa frequência dominam; aparência de "memória persistente".


# Apêndice B — Formulário essencial

Tabela de bolso com as fórmulas que aparecem ao longo do documento. Contexto e hipóteses são descritos nas seções correspondentes; aqui o objetivo é consulta rápida.

### B.1 Distribuições fundamentais

| Quantidade | Fórmula | Onde aparece |
|---|---|---|
| Binomial $\mathrm{Bin}(n,p)$ | $P(X=k)=\binom{n}{k}p^k(1-p)^{n-k}$ | contagens de sucessos em $n$ ensaios |
| Poisson $\mathrm{Pois}(\mu)$ | $P(X=k)=e^{-\mu}\mu^k/k!$; $E[X]=\mu$, $\mathrm{Var}(X)=\mu$ | eventos raros por unidade de tempo/área |
| Geométrica $\mathrm{Geom}(p)$ | $P(X=k)=(1-p)^{k-1}p$; $E[X]=1/p$ | espera até o primeiro sucesso |
| Exponencial $\mathrm{Exp}(\lambda)$ | densidade $\lambda e^{-\lambda x}$; sem memória: $P(X>x+a \mid X>a)=P(X>x)$ | inter-chegadas do Poisson; tempos de serviço |
| Normal padrão | $\Phi(z)=P(Z\le z)$; tabela usual | limites, intervalos de confiança |
| Uniforme $\mathrm{Unif}(a,b)$ | $E[X]=(a+b)/2$, $\mathrm{Var}=(b-a)^2/12$ | erros de arredondamento, prioris |

Somas úteis: se $X_1,\ldots,X_n$ são independentes, $E[\sum X_i]=\sum E[X_i]$ e $\mathrm{Var}(\sum X_i)=\sum \mathrm{Var}(X_i)$; se todas têm a mesma distribuição com média $\mu$ e variância $\sigma^2$, então $\bar X_n$ tem variância $\sigma^2/n$ — a raiz quadrada disso é o erro-padrão $\sigma/\sqrt n$.

### B.2 Passeio aleatório simples

Definição: $S_n = Z_1+\cdots+Z_n$ com $P(Z_i=\pm1)=p,q=1-p$ i.i.d.

- Média e variância: $E[S_n]=n(p-q)=n(2p-1)$; $\mathrm{Var}(S_n)=4npq$ (no caso simétrico, $\mathrm{Var}=n$).
- Recorrência/transiência (Pólya): passeio **simétrico** em redes inteiras $\mathbb{Z}^d$ volta à origem com probabilidade 1 se $d\le 2$; se $d\ge 3$, a probabilidade é menor que 1 (transiente).
- Ruína do jogador (com dois patrimônios finitos, jogo justo): $P(\text{atingir } N \text{ antes de } 0 \mid \text{começa em } i) = i/N$. Com viés $p\ne q$: $P_{\text{ruina}} = \dfrac{(q/p)^i-(q/p)^N}{1-(q/p)^N}$.
- Escala (intuição do movimento browniano): $W^{(n)}_t=\frac{1}{\sqrt{n}}S_{\lfloor nt\rfloor}$ converge em lei para um movimento browniano quando $n\to\infty$ (Donsker).

### B.3 Cadeias de Markov

- Discretas, matriz $P$ (linha estocástica): $P_{ij}=P(X_{n+1}=j \mid X_n=i)$; Chapman-Kolmogorov: $P^n$ dá probabilidades de $n$ passos.
- Distribuição estacionária: $\boldsymbol{\pi}P=\boldsymbol{\pi}$, $\sum_i \pi_i=1$; solução única para cadeias irredutíveis e aperiódicas (regulares); convergência $\lim_{n\to\infty}P^n$ existe nessas condições.
- Estados absorventes: $i$ absorvente se $P_{ii}=1$; probabilidades de absorção resolvem sistema linear harmônico $(I-P_Q)p=h$ com condições de contorno adequadas.
- Contínuas (CTMC): matriz de taxas $Q$ (linhas somando 0); equações de Kolmogorov: $P'(t)=P(t)Q$ (direta) e $P'(t)=QP(t)$ (inversa); $\boldsymbol{\pi}Q=0$ para distribuição estacionária.

### B.4 Processo de Poisson

- Definição (via incrementos): $N(0)=0$; incrementos independentes; incrementos estacionários; $N(t)\sim\mathrm{Pois}(\lambda t)$.
- Inter-chegadas: $T_k-T_{k-1}$ i.i.d. $\mathrm{Exp}(\lambda)$; $E[T_k]=1/\lambda$.
- Propriedades: $E[N(t)]=\mathrm{Var}(N(t))=\lambda t$; superposição (processos independentes somam taxas $\lambda_1+\lambda_2$); afinamento (manter cada evento com prob. $p$ mantém forma Poisson com taxa $\lambda p$).
- Generalizações: Poisson não-homogêneo ($\lambda(t)$ determinístico), Poisson composto, campos pontuais de Poisson espaciais.

### B.5 Movimento browniano e cálculo estocástico

- Padrão: $W_0=0$; incrementos independentes; $W_t-W_s\sim N(0,t-s)$; caminhos contínuos; $E[W_t]=0$, $\mathrm{Var}(W_t)=t$; covariância $\mathrm{Cov}(W_s,W_t)=\min(s,t)$.
- Martingalas associadas: $W_t$ e $W_t^2-t$ são martingalas; variação quadrática $[W]_t=t$.
- Versão com deriva/volatilidade: $X_t=\mu t+\sigma W_t$; Geometric Brownian motion: $S_t=S_0\exp((\mu-\sigma^2/2)t+\sigma W_t)$.
- Precificação (medida neutra ao risco): $C_0=e^{-rT}E[(S_T-K)^+]$; Black-Scholes: $C=S\Phi(d_1)-Ke^{-rT}\Phi(d_2)$, $d_{1,2}=\frac{\ln(S/K)+rT}{\sigma\sqrt T}\mp\frac{\sigma\sqrt T}{2}$.

### B.6 Renovação

- Processo de renovação: instantes $T_k=\tau_1+\cdots+\tau_k$ com $\tau_i$ i.i.d. comuns; $N(t)=\max\{k:T_k\le t\}$.
- Teorema elementar da renovação: $E[N(t)]\sim t/E[\tau_1]$ (para $F$ não-aritmética; notação assintótica usual).
- Paradoxo da inspeção: intervalo observável em um instante escolhido aleatoriamente é maior, em média: $E[\text{intervalo observado}]=E[\tau_1^2]/E[\tau_1]\ge E[\tau_1]$.
- Vida residual contínua: $P(R>l+s\mid R>l)=\dfrac{1-F(l+s)}{1-F(l)}$.
- Recompensa renovação: taxa média $= \dfrac{E[\text{recompensa por ciclo}]}{E[\text{duração do ciclo}]}$ (base de modelos econômicos de manutenção e seguros).

### B.7 Séries temporais

- Estacionariedade estrita: distribuição conjunta invariante por translação temporal; fraca: $E[X_t]$ constante e $\mathrm{Cov}(X_s,X_t)=\gamma(s-t)$.
- FAC: $\rho(k)=\gamma(k)/\gamma(0)$; propriedades: positiva definida, definida positiva; AR(1) cai exponencialmente; MA(j) corta após $j$ lags.
- ARIMA(p,d,q): modelo $(1-\sum_{i=1}^p\phi_i B^i)(1-B)^d X_t = (1+\sum_{j=1}^q\theta_j B^j)\varepsilon_t$ com ruído branco $\varepsilon_t$; condições de invertibilidade e causalidade pelas raízes do polinômio.
- Diagnóstico: Ljung–Box testa ausência de autocorrelação nos resíduos; testes ADF/KPSS para raízes unitárias; transformações para estabilizar variância (log, Box-Cox).
- Previsão: intervalos preditivos ampliam com o horizonte; backtesting obrigatório.

### B.8 Filas (M/M/1) e epidemias (SIR)

M/M/1 com chegada Poisson($\lambda$), serviço Exp($\mu$), $\rho=\lambda/\mu<1$:

$$P_n=(1-\rho)\rho^n,\quad L_q=\frac{\rho^2}{1-\rho},\quad W_q=\frac{\rho^2}{\mu(1-\rho)},\quad W=W_q+\frac1\mu,\quad L=\lambda W\ \ (\text{lei de Little}).$$

SIR: $\dot S=-\beta SI/N$, $\dot I=\beta SI/N-\gamma I$, $\dot R=\gamma I$; número básico $R_0=\beta N/\gamma$ (taxas per capita constantes, população homogênea). Limiar: $R_0<1$ extingue, $R_0>1$ gera surto (probabilisticamente na versão estocástica). Para populações grandes com suscetível inicial próximo de 1, tamanho final aproximado pela relação implícita clássica $\ln s_\infty\approx-R_0(1-s_\infty)$ (com fração infectada inicial pequena); resolver iterativamente.

### B.9 Simulação e inferência

- Erro Monte Carlo: estimador de média com $N$ replicações independentes tem desvio típico $\sigma/\sqrt N$ — quadruplicar amostras melhora precisão por fator 2.
- Redução de variância: antitéticos ($\mathrm{Var}(\bar Y)\approx \mathrm{Var}(Y)/2$ quando a correlação entre pares é próxima de $-1$), controles variáveis (usar correlação com variável conhecida), estratificação.
- p-valor: $p=P(\text{estatística} \ge \text{observada} \mid H_0)$; interpretação correta = medida de incompatibilidade com $H_0$ sob o modelo; **não** é $P(H_0)$ nem tamanho do efeito.
- Correções de multiplicidade: Bonferroni ($\alpha/m$ para $m$ testes), FDR (Benjamini-Hochberg), permutação/min-P; pré-registro resolve problemas estruturais que correções post-hoc não resolvem.
- IC 95% frequente: cobertura a longo prazo de repetições hipotéticas — cuidado ao interpretar um único intervalo; alternativas: bootstrap, Bayesian credible intervals.

---


---

# Apêndice C — Referências citadas

Compilação automática de todas as fontes com hyperlink usadas neste documento,
agrupadas por tipo. Consulte também o Capítulo 8 para trilhas de estudo.

### Materiais brasileiros (universidades)

- [Curso IPRJ/PUC-Rio](https://www.maxwell.vrac.puc-rio.br/15321/15321_4.PDF)
- [ME501 — Unicamp](https://www.ime.unicamp.br/graduacao/disciplina/me501-processos-estocasticos)
- [Notas de processos estocásticos IME-USP](https://www.ime.usp.br/~lane/home/MAE5870/STprocestv2.pdf)
- [slides UFSC](https://www.inf.ufsc.br/~marcelo.menezes.reis/Processos01.pdf)
- [Reprodutibilidade na ciência — USP](https://cienciaaberta.usp.br/reprodutibilidade-na-ciencia/)
- [Página oficial](https://www.ime.unicamp.br/pos-graduacao/disciplina/mi414-introducao-aos-processos-estocasticos)
- [Registro JupiterWeb](https://uspdigital.usp.br/jupiterweb/obterDisciplina?sgldis=mae0312)
- [PDF](https://www.ime.usp.br/~galves/cursos/NotasdeAulaAplProcEst.pdf)
- [PDF](https://www.est.ufmg.br/portal/wp-content/uploads/2023/01/RTE-01-2011.pdf)
- [PDF](https://www.professor.ufop.br/sites/default/files/magno/files/capitulo_6_-_processos_estocasticos_e_cadeias_de_markov_1.pdf)
- [PDF](https://engcomp.ufc.br/wp-content/uploads/2020/08/ti0112-processos-estocasticos.pdf)
- [Emap](https://emap.fgv.br/disciplina/processos-estocasticos-1)
- [See](https://see.cepel.br/manual/anatem/recursos/serie/estocastico.html)

### Enciclopédias e verbetes

- [Stochastic process](https://en.wikipedia.org/wiki/Stochastic_process)
- [0, \infty)$).
- $S$: onde o processo vive. Pode ser finito (ex.: "saudável/doente"), enumerável (ex.: contagens $\{0, 1, 2, \dots\}$) ou contínuo (ex.: posições na reta).

Note que a palavra "processo" **não exige** que o índice seja temporal: um campo aleatório em uma imagem é indexado por coordenadas de pixels, e uma série de Fourier aleatória por frequências. Mas o uso esmagadoramente dominante é como modelo de evolução temporal.

### 1.3 Três exemplos construtivos (para fixar a ideia)

**(a) Sequência de moedas honestas.** Jogue uma moeda justa repetidamente. O espaço de resultados $\Omega$ pode ser pensado como o conjunto de sequências infinitas de caras ($H$) e coroas ($T$). Defina $X_n = 1$ se a $n$-ésima jogada for cara e $0$ caso contrário. Cada $X_n$ é uma cópia de um Bernoulli$(1/2)$, e juntos formam um processo estocástico em tempo discreto com espaço de estados $\{0, 1\}$. Este é o "átomo" de quase todos os processos: muitos modelos são sequências de decisões tipo cara/coroa com regras de transição mais elaboradas.

**(b) Posição de uma partícula na reta (passeio aleatório).** Uma partícula começa em $0$ e, a cada passo, anda $+1$ ou $-1$ com igual probabilidade, independentemente do passado. A posição no passo $n$ é $S_n = Z_1 + \cdots + Z_n$ com $Z_i \in \{-1, +1\}$. Aqui $T=\{0,1,2,\dots\}$ e $S=\mathbb{Z}$ (espaço de estados infinito enumerável). Esse processo será nosso primeiro grande personagem na seção 4.1.

**(c) Contagem de eventos até o instante $t$ (processo de contagem).** Considere chamadas telefônicas chegando a uma central. Para cada $t \geq 0$, defina $N(t)$ = "número de chamadas chegadas entre $0$ e $t$". Cada $N(t)$ é uma variável de contagem; o vetor inteiro $(N(t))_{t \geq 0}$ é um processo com tempo contínuo e espaço de estados $\{0, 1, 2, \dots\}$. Quando impomos axiomas razoáveis de independência e taxa constante, obtemos o **processo de Poisson** (seção 4.3).

### 1.4 Duas leituras que você precisa separar

Esta é a distinção conceitual mais importante da seção — e fonte de muita confusão:

1. **Instante fixo, resultado variável:** fixe $t$ e olhe $X_t$ como variável aleatória. Perguntamos: "qual a distribuição do estado no instante $t$?" Ex.: qual a chance de chover dia 10? (Ou: de haver 3 clientes na fila às 15h?)
2. **Resultado fixo, tempo variável:** fixe uma realização $\omega$ (o "sorteio" específico do universo) e olhe a função $t \mapsto X_t(\omega)$. Essa função é a **trajetória**, **caminho amostral** ou **realização** do processo. Ex.: a série histórica real de chuva do ano passado é uma trajetória.

Um gráfico de uma simulação mostra **uma trajetória**; o processo em si é o *mecanismo probabilístico* que gera todas as trajetórias possíveis. Você nunca vê o processo inteiro diretamente — só amostras dele. Como resumem as notas da UFABC: "um processo estocástico pode ser pensado como uma regra aleatória que escolhe uma trajetória; antes de sabermos qual ocorreu, temos várias trajetórias possíveis; depois que $\omega$ é fixado, uma delas fica determinada".

Quando quiser comunicar isso em uma frase: *a variável aleatória responde "onde podemos estar"; a trajetória conta "como foi desta vez"*.

### 1.5 De onde vêm os processos? A construção de Kolmogorov

Como se *constrói* um processo de fato? Uma abordagem moderna (devida a Andrey Kolmogorov) diz: especifique primeiro, para cada número finito de instantes $t_1, \dots, t_n$, a distribuição conjunta de $(X_{t_1}, \dots, X_{t_n})$ — as chamadas **distribuições de dimensão finita** —, e faça isso de forma **consistente**: a distribuição que você especifica para $n+1$ instantes deve ter como marginais as distribuições que você já tinha para $n$ instantes (a "condição de consistência"). Sob condições adequadas, existe então uma única lei de probabilidade que realiza todas essas distribuições simultaneamente — isto é, um processo com essas leis de dimensão finita. Esse é o conteúdo do [teorema de extensão de Kolmogorov](https://en.wikipedia.org/wiki/Kolmogorov_extension_theorem)
- [Central limit theorem](https://en.wikipedia.org/wiki/Central_limit_theorem)
- [Gambler's ruin — Wikipédia](https://en.wikipedia.org/wiki/Gambler%27s_ruin)
- [Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process)
- [verbete Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
- [Renewal theory](https://en.wikipedia.org/wiki/Renewal_theory)
- [compartmental models](https://en.wikipedia.org/wiki/Compartmental_models_(epidemiology)
- [M/M/1 queue — Wikipedia](https://en.wikipedia.org/wiki/M/M/1_queue)

### Notas e cursos internacionais

- [An Introduction to Stochastic Processes in Continuous Time](https://www.math.leidenuniv.nl/~spieksma/colleges/sp-master/LN.pdf)
- [Law of large numbers and the Central Limit Theorem — notas do MIT](https://math.mit.edu/~dav/05.dir/class6-prep.pdf)
- [notes de aula CMU, "Random Times and Their Properties"](https://www.stat.cmu.edu/~cshalizi/754/notes/lecture-06.pdf)
- [resumo didático](https://math.dartmouth.edu/~pw/math100w13/mare.pdf)
- [MIT 6.262, capítulo 2](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/3a19ce0e02d0008877351bfa24f3716a_MIT6_262S11_chap02.pdf)
- [cap. 4 — Renovação](https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/931ffa0940899c27f34b71ad64fd2bb0_MIT6_262S11_chap04.pdf)
- [lecture 2](https://www.stat.cmu.edu/~cshalizi/754/notes/lecture-02.pdf)
- [BM](https://galton.uchicago.edu/~lalley/Courses/313/BrownianMotionCurrent.pdf)
- [martingalas](http://galton.uchicago.edu/~lalley/Courses/385/Martingales.pdf)
- [PDF](https://www.math.purdue.edu/~stindel/teaching/ma539/brownian-motion1.pdf)
- [teorema da parada opcional (LaLonde)](https://math.dartmouth.edu/~pw/math100w13/lalonde.pdf)

### R: pacotes, documentos e comunidades

- [CRAN Task View: Time Series Analysis](https://cran.r-project.org/view=TimeSeries)
- [Getting started with spatstat](https://cran.r-project.org/web/packages/spatstat/vignettes/getstart.pdf)
- [CRAN SpatioTemporal task view](https://cran.r-project.org/view=SpatioTemporal)
- [vinheta oficial RStan](https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html)
- [Slides](https://curso-r.github.io/main-intro-programacao/slides/01-introducao-ao-curso.pdf)
- [Repositório](https://github.com/arthurspk/guiader)
- [lista completa](https://cran.r-project.org/web/views/)
- [Finance](https://cran.r-project.org/view=Finance)

### Metodologia, p-valores e reprodução científica

- [The Statistical Crisis in Science — American Scientist](https://www.americanscientist.org/article/the-statistical-crisis-in-science)
- [Big little lies: a compendium and simulation of p-hacking — Royal Society Open Science](https://royalsocietypublishing.org/rsos/article/10/2/220346/92017/Big-little-lies-a-compendium-and-simulation-of-p)
- [Impact of redefining statistical significance on p-hacking — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11098386/)
- [Error bars in experimental biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC2064100/)
- [Robust misinterpretation of confidence intervals — PMC](https://link.springer.com/article/10.3758/s13423-013-0572-3)
- [On the Assessment of Monte Carlo Error — Am. Statistician](https://pmc.ncbi.nlm.nih.gov/articles/PMC3337209/)

### Artigos e revistas brasileiras

- [Lau et al., TEMA/SBMAC](https://tema.sbmac.org.br/tema/article/download/1670/1173)
- [PDF](https://www.sbmac.org.br/wp-content/uploads/2022/08/livro_81.pdf)
- [SciELO](https://www.scielo.br/j/rac/a/FyNM3rThZN7wQgMKshwLHjk/?lang=en)

### Blogs, fóruns e tutoriais comunitários

- [por que em 1D/2D a recusa é certa](https://math.stackexchange.com/questions/3233552/intuitive-reason-why-a-simple-symmetric-random-walk-is-recurrent-on-bbb-z2-a)
- [discussão StackExchange sobre a variação quadrática do Wiener process](https://math.stackexchange.com/questions/3705345/computation-of-the-quadratic-variation-of-wiener-process)
- [Random Seeds and Reproducibility](https://medium.com/data-science/random-seeds-and-reproducibility-933da79446e3)
- [Choosing the "Correct" Seed for Reproducible Research, StackExchange](https://stats.stackexchange.com/questions/335936/choosing-the-correct-seed-for-reproducible-research-results)
- [Stack Overflow: reasons for using set.seed](https://stackoverflow.com/questions/13605271/reasons-for-using-the-set-seed-function)
- [Setting a seed in R when using parallel simulation](https://www.r-bloggers.com/2018/07/%F0%9F%8C%B1-setting-a-seed-in-r-when-using-parallel-simulation/)
- [exemplo implementado em R, StackExchange Quant](https://quant.stackexchange.com/questions/50168/monte-carlo-option-pricing-with-r)
- [Texto](https://medium.com/@flaviagaia/introdu%C3%A7%C3%A3o-%C3%A0s-s%C3%A9ries-temporais-com-a-linguagem-r-37157435cfbd)
- [Publicação](https://rpubs.com/RobsonOL/series_temporais_r)

### Outros recursos

- [Introdução aos Processos Estocásticos (UFABC)](https://danielmiranda.prof.ufabc.edu.br/livros/processos/capitulo-01-processos-estocasticos.html)
- [discussão no Risk Hub](https://riskhub.org/blogs/kolmogorovs-extension-perspective-on-stochastic-processes)
- [basic concepts of the Poisson process](https://www.probabilitycourse.com/chapter11/11_1_2_basic_concepts_of_the_poisson_process.php)
- [filtrations and stopping times](https://www.randomservices.org/random/prob/Stop.html)
- [Martingales and the Optional Stopping Theorem](https://www.jeremykun.com/2014/03/03/martingales-and-the-optional-stopping-theorem/)
- [notes of Marc Yor / Samuel Drapeau](https://www.samuel-drapeau.info/downloads/WS2015/StochProcess/Stochastic_Processes_chap02_1126.pdf)
- [notas introdutórias](https://www.cmat.edu.uy/~lessa/resource/randomwalknotes.pdf)
- [Por que passeios aleatórios se perdem em 3D](https://www.youtube.com/watch?v=byvEzyFgv44)
- [nota histórica](http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-GR.pdf)
- [derivação via funções geradoras, Grinstead & Snell](https://stats.libretexts.org/Bookshelves/Probability_Theory/Introductory_Probability_(Grinstead_and_Snell)
- [tutorial espectral para cientistas de dados](https://arxiv.org/pdf/2207.02296)
- [estações estacionárias e limites](https://www.probabilitycourse.com/chapter11/11_3_2_stationary_and_limiting_distributions.php)
- [nota introdutória sobre variação quadrática em finanças](https://www.finance-tutoring.fr/the-relation-between-brownian-motion-and-quadratic-variation-in-layman%E2%80%99s-terms%E2%80%A6?mobile=1)
- [notas de renovação da Columbia](http://www.columbia.edu/~ks20/4106-18-Fall/Notes-Renewal-Theory.pdf)
- [MathWorks — Time Series Regression IV: Spurious Regression](https://www.mathworks.com/help/econ/time-series-regression-iv-spurious-regression.html)
- [estudo sobre não-estacionaridade](https://mpra.ub.uni-muenchen.de/27926/1/Stationarity_of_time_series_and_the_problem_of_spurious_)
- [notas CGSpace/ CGIAR sobre pitfalls de raiz unitária](https://cgspace.cgiar.org/server/api/core/bitstreams/7aa92438-0578-4751-b149-6a8f8f12dfb8/content)
- [Spatial Statistics for Data Science](https://www.paulamoraga.com/book-spatial/the-spatstat-package.html)
- [tutorial gstat/INPE](http://www.dpi.inpe.br/gilberto/tutorials/software/R-contrib/gstat.pdf)
- [capítulo Kriging do livro de Moraga](https://www.paulamoraga.com/book-spatial/kriging.html)
- [documentação de MCMC do Stan](https://mc-stan.org/docs/reference-manual/mcmc.html)
- [Monte Carlo Convergence and Standard Error, Risk Hub](https://riskhub.org/blogs/monte-carlo-convergence-and-standard-error)
- [Variance Reduction Techniques for Efficient Monte Carlo](https://mbrenndoerfer.com/writing/variance-reduction-monte-carlo-simulation)
- [Guias práticos para uma pesquisa mais reprodutível](https://www.reprodutibilidade.org/post/guias-praticos-para-uma-pesquisa-mais-reprodutivel)
- [Black–Scholes vs Monte Carlo na prática](https://www.acturtle.com/blog/call-option-valuation-black-scholes-vs-monte-carlo)
- [tutorial do King sobre simulação estocástica, com código R no pacote pomp](https://kingaa.github.io/short-course/stochsim/stochsim.html)
- [guia aplicado M/M/1 vs M/M/c](https://metricgate.com/blogs/queueing-theory-mm1-mmc-models/)
- [texto completo](https://sites.stat.columbia.edu/gelman/research/unpublished/p_hacking.pdf)
- [Time Series Nonstationarity — bookdown R for Political Scientists](https://bookdown.org/ccolonescu/RPoE4/time-series-nonstationarity.html)
- [INFORMS/WSC 2001 paper](https://informs-sim.org/wsc01papers/006.PDF)
- [Monte Carlo Simulation — practical guide](https://mcpanalytics.ai/articles/monte-carlo-simulation-practical-guide-for-data-driven-decisions)
- [Fama — Random Walks in Stock-Market Prices](https://www.chicagobooth.edu/~/media/34F68FFD9CC04EF1A76901F6C61C0A76.PDF)
- [Investopedia — Random Walk Theory](https://www.investopedia.com/terms/r/randomwalktheory.asp)
- [Samuelson vs Fama on the Efficient Market Hypothesis — Œconomia](https://journals.openedition.org/oeconomia/5300)
- [Página do curso](https://www.escolavirtual.gov.br/curso/325)
- [Alura](https://www.alura.com.br/artigos/linguagem-r)
- [Análise Macro](https://analisemacro.com.br/curso/analise-de-series-temporais-usando-r/)
- [epiR Handbook (versão PT)](https://www.epirhandbook.com/pt/new_pages/time_series.pt.html)
- [introdução à renovação](https://www.randomservices.org/random/renewal/Introduction.html)

# Apêndice D — Soluções comentadas (esboços)

**Seção 1.**
1. Uma trajetória fixa mostra apenas **um** desfecho possível; a distribuição em um instante fixo captura a incerteza naquele momento. Duas leituras complementares, como na definição formal.
2. Ajustar uma curva a um caminho simulado tende a overfitting do ruído: você estaria modelando uma realização específica, não a lei geradora. Valide com múltiplas trajetórias ou dados out-of-sample.

**Seção 2.**
1. $X_i\sim\mathrm{Bernoulli}(p)$: esperança $p$, variância $p(1-p)$. Como soma de $n$ idênticas independentes, $E[\sum X_i]=np$, $\mathrm{Var}(\sum X_i)=np(1-p)$ — recupera binomial como caso concreto.
2. Na decomposição do problema da moeda honesta, a variância cai linearmente com $n$ porque a dependência entre os fatores de escala não altera a independência das variáveis subjacentes. Verifique com simulações diretas.
3. Intuição via TCL: qualquer agregado com muitas componentes pequenas e fracamente dependentes se aproxima de uma Normal — daí sua ubiquidade.

**Seção 3.**
1. Um tempo de parada exige decidir "parar agora?" usando somente informação até hoje, sem prever o futuro; o máximo de um processo até um alvo é exemplo clássico.
2. Jogo justo: $E[X_{n+1}\mid \mathcal{F}_n]=X_n$; o capital acumulado de um apostador de sorte justa é martingala; para jogos com vantagem do cassino, é supermartingala (esperança decrece).
3. Condição opcional do teorema exige integrabilidade uniforme — falha quando o jogo pode prolongar-se indefinidamente com apostas que dobram (paradoxos tipo doubling); sempre verificar.

**Seção 4.1 (Passeio aleatório).**
1. Simulação: conte retornos à origem para 1D/2D/3D; nas duas primeiras dimensões o número médio cresce muito (recorrente), na terceira fica limitado (transiente, retorno < 1). Compare com a teoria de Pólya.
2. $E[S_n]=0$, $\mathrm{Var}(S_n)=n$ no caso simétrico (cada passo tem variância 1). Use para conferir as faixas típicas $\pm\sqrt n$ dos gráficos.
3. O passeio assimétrico $p>q$ tem deriva $2p-1>0$: vai a $+\infty$ com probabilidade 1, logo não é recorrente mesmo em 1D. Isso generaliza: recorrência exige balanço detalhado, não só simetria de estados.
4. Ruína do jogador: monte as equações $P_i=(P_{i-1}+P_{i+1})/2$ com condições de contorno $P_0=0,P_N=1$; a solução é $P_i=i/N$. Com viés, derive a fórmula geral via razões $(q/p)^i$ e discuta limite $N\to\infty$ (probabilidade de ruína 1 para qualquer drift adverso).

**Seção 4.2 (Cadeias de Markov).**
1. Resolva $(I-P^\top)\boldsymbol{\pi}=\boldsymbol{0}$ com restrição $\sum\pi_i=1$; numérico: autovetor de $P^\top$ com autovalor mais próximo de 1. Confirme $\boldsymbol{\pi}P=\boldsymbol{\pi}$ e que frequências simuladas longas se aproximam de $\boldsymbol{\pi}$.
2. Se todos os estados se comunicam e pelo menos um retorna em tempos com GCD 1, então cadeia regular; teste elevando $P$ a potências crescentes e observando convergência das linhas.
3. Ajuste de modelo: escolha entre 2 ou 3 estados segundo critérios interpretáveis (tempo de residência) e valide com testes de razão de verossimilhança/BIC em dados.

**Seção 4.3 (Poisson).**
1. Gere `chegadas <- cumsum(rexp(...))` e compare histogramas de $N(t)$ empírico com pmf teórica; confirme a memória praticamente nula (correlação série de contagens em intervalos distintos quase nula).
2. Superposição de dois Poissons com taxas 2 e 3 dá Poisson(5) — bom teste de sanidade numérico.
3. Teste de qualidade: gráfico quantílico (QQ) exponencial das inter-chegadas versus linha teórica; desvios sistemáticos indicam renovação (distribuição não-exponencial).

**Seção 4.4 (Browniano).**
1. Empírico: rode 200 simulações com passos pequenos, calcule a variância em tempos fixos e ajuste contra $t$ — inclinação deve estar próxima de 1.
2. Construa $X_t=\mu t+\sigma W_t$ por composição direta; confira média e variância empíricas.
3. Compare caminhos do passeio escalado $W^{(n)}$ com BM: quanto menor o passo de tempo, mais próximo visualmente; a diferença concentra-se em saltos grandes (que o passeio não tem).
4. Reflexão: para o maximum de BM até $t$, vale $P(\sup_{s\le t}W_s\ge a)=2P(W_t\ge a)$ — princípio da reflexão; útil para derivar distribuições de máximos em barreiras.

**Seção 4.5 (Renovação).**
1. Simule inter-chegadas $ \mathrm{Exp}(1)$ (média 1) e $\mathrm{Gama}(2,1)$ (média 2); conte $N(t)$ até $t=2000$; compare $E[N(t)]/t$ com $1/E[\tau]$ obtendo algo próximo de 1 e 0,5 respectivamente.
2. Demonstre o paradoxo medindo intervalos observados num instante fixo versus intervalos típicos das renovais; o primeiro será substancialmente maior.
3. Modelagem de confiabilidade: se componente envelhece (failure rate crescente), Poisson subestima falhas futuras; use renovação com Weibull/lognormal.

**Seção 4.6 (Séries temporais).**
1. Sobre a chuva anual paulistana (dados reais da Seção 6.D): lag-1 ≈ 0,36, φ OLS ≈ 0,34; teste ADF deve rejeitar raiz unitária forte (confirmar), e diferenciação produz série branca com FAC sem estrutura. Decadais mostram leve tendência negativa — exija testes formais antes de afirmar.
2. Para mensais, SARIMA$(0,1,0)(0,1,1)_{12}$ combina diferenciação ordinária (nível) e sazonal; cheque resíduos com Ljung-Box e faça backtest curto.
3. Long memory: se Hurst > 0,5 em análise R/S, considere modelagem fracionária; cuidado com amostras curtas.

**Seção 5.**
1. Repita o experimento antitético para $N=10^4$ pares vs. métodos comuns; a queda de variância deve ser visível (ordem de metade, dependendo do estimando).
2. Documente versões (`sessionInfo()`), ambiente (`renv::snapshot()`) e sementes; registre todas as sementes usadas.
3. Faça um teste com três códigos de discretização de BM (Euler vs. Milstein vs. salto composto) e discuta custos/tempo.

**Seção 6.**
1. Opções: com $\sigma=0{,}4$, o erro MC deve subir substancialmente (payoff dispersa mais) e a convergência relativa piora; compare com BS para confirmar. Sanity checks devem passar.
2. SIR: plote $\hat p_{\text{extinção}}(R_0)$ com 500 simulações por ponto; a curva cai rapidamente perto de $R_0=1$; para $R_0\ll 1$, extinção quase certa.
3. M/M/1 ($\lambda=4,\mu=6$): $\rho=2/3$; $W_q=(4/9)/(6\cdot1/3)=2/9\approx0{,}222$ min; $L_q=\lambda W_q=8/9\approx0{,}889$; $W=W_q+1/\mu=17/36\approx0{,}472$ min; $L=\lambda W=17/9\approx1{,}889$ — identidade $L=\lambda W$ exata na teoria; simulação confirma dentro do erro Monte Carlo.
4. Clima: a FAC das primeiras diferenças deve revelar poucos lags significativos, mas checar quadrantes/meses (a sazonalidade é forte na série mensal; diferenciar anualmente não remove a estrutura intra-anual).

**Seção 7.**
1. Com oito testes, $\Pr(\min p<0{,}05)$ sobe ainda mais acima de 0,5; Bonferroni usaria $0{,}05/8=0{,}00625$ — conservador porque assume independência; a estrutura real dos testes reduz o custo efetivo (multiverse/permutation calibration).
2. Regressão entre passeios independentes: o R² tipicamente fica alto (tendências compartilhadas) mesmo com coeficiente "significativo"; séries relacionadas produzem R² semelhante — ilustra regressão espúria. Sempre teste estacionariedade antes de regredir.
3. AR(1), $\phi=0{,}95$, $n=50$: OLS $\hat\phi$ enviesado para baixo (para 0); a inferência baseada em $|\hat\phi|>0$ seria enganosa sobre a força real. Amostragem bayesiana/bootstrap corrige parcialmente; reporte intervalos.
4. Bars de erro de média ($\sigma/\sqrt n$) ficam mais estreitas que intervalos de predição porque só descrevem a incerteza da **média**, não de uma nova observação; use barras de previsão quando quiser cobrir valores individuais futuros; nunca trate "barras que se tocam" como prova de nulidade (Cumming).

---

# Sobre este documento

Esta documentação foi compilada em novembro de 2026 como material de estudo sobre processos estocásticos, com ênfase em intuição, rigor moderado e implementação em **R**. As fontes citadas (notas universitárias, livros abertos, artigos, páginas oficiais de cursos e bases de dados públicas consultadas diretamente) estão listadas no Apêndice C, junto com sugestões de trilha de estudo no Capítulo 8. **Aviso:** os scripts R são ferramentas didáticas; resultados dependem de sementes, parâmetros e dados — sempre reproduza localmente e valide contra seus próprios dados antes de tomar decisões.

*Sugestões de melhoria (correções, exemplos novos, traduções de notação) são bem-vindas — mantenha este documento vivo como parte do seu caderno de estudos.*


---

