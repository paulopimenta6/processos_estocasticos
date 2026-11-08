# =============================================================================
# demo_processos_estocasticos.R
# Script de apoio da documentação "Processos Estocásticos: guia ilustrado com R"
# =============================================================================
# Uso: rode em R >= 4.x. Blocos são didáticos e independentes entre si.
# Pacotes opcionais (instale conforme a seção que deseja reproduzir):
#   install.packages(c("gstat", "spatstat", "rstan"))   # seções 5.4-5.5
# Para séries temporais e clima usamos apenas pacotes básicos (stats, graphics).
# Atenção: resultados dependem das sementes fixadas em cada bloco (set.seed).
# =============================================================================


# -------------------------------------------------------------------------- #
# Bloco 1 (extraído da documentação)
# -------------------------------------------------------------------------- #
set.seed(42)
n <- 200                 # tamanho de cada amostra
nsim <- 1000             # número de repetições do experimento
medias <- replicate(nsim, mean(sample(1:6, size = n, replace = TRUE)))
mean(medias)             # deve estar perto de 3.5
sd(medias)               # ≈ sqrt(35/12)/sqrt(n) ≈ 0.17

# -------------------------------------------------------------------------- #
# Bloco 2 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 3 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 4 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 5 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 6 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 7 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 8 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 9 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 10 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 11 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 12 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 13 (extraído da documentação)
# -------------------------------------------------------------------------- #
set.seed(2026)             # fixa a semente da sequência inteira
amostra <- sample(1:6, size = 100, replace = TRUE)   # 100 lançamentos de dado
media_sim <- mean(amostra)

# Repita a linha acima quantas vezes quiser: mesmo resultado sempre.

# -------------------------------------------------------------------------- #
# Bloco 14 (extraído da documentação)
# -------------------------------------------------------------------------- #
# Média empírica da soma de 100 dados — converge para 350 pela LGN
sim <- replicate(10_000, sum(sample(1:6, size = 100, replace = TRUE)))
hist(sim, breaks = 30, freq = FALSE, main = "Soma de 100 dados: 10.000 simulações")
abline(v = 350, col = "red", lwd = 2)
c(media = mean(sim), desvio_padrao = sd(sim))

# -------------------------------------------------------------------------- #
# Bloco 15 (extraído da documentação)
# -------------------------------------------------------------------------- #
set.seed(7)
serie <- arima.sim(list(order = c(2, 1, 1), ar = c(0.8, -0.2), ma = 0.3), n = 500)
plot(serie, main = "Série sintética ARIMA(2,1,1) — 500 pontos")

# -------------------------------------------------------------------------- #
# Bloco 16 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 17 (extraído da documentação)
# -------------------------------------------------------------------------- #
# if (!requireNamespace("spatstat", quietly = TRUE)) install.packages("spatstat")
library(spatstat)
win <- owin(xrange = c(0, 1), yrange = c(0, 2))       # janela de observação
X <- runifpoint(n = 200, win = win)                   # 200 pontos uniformes
plot(X, main = "Padrão pontual simulado (uniforme)")
intensity(X)                 # intensidade média estimada
nn <- nndist(X)              # distâncias ao vizinho mais próximo
hist(nn$dist, main = "Distância ao vizinho mais próximo")

# -------------------------------------------------------------------------- #
# Bloco 18 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 19 (extraído da documentação)
# -------------------------------------------------------------------------- #
# if (!requireNamespace("rstan", quietly = TRUE)) install.packages("rstan")
# library(rstan)
# fit <- stan(file = "modelo.stan", data = list(N = length(y), y = y),
#             iter = 2000, warmup = 1000, chains = 4, seed = 2026)
# print(fit)                      # resumo posterior: média, dp, ICs 95%
# summary(fit)$summary[, "Rhat"]  # diagnóstico de convergência (deve ser ~1)

# -------------------------------------------------------------------------- #
# Bloco 20 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 21 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 22 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 23 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 24 (extraído da documentação)
# -------------------------------------------------------------------------- #
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

# -------------------------------------------------------------------------- #
# Bloco 25 (extraído da documentação)
# -------------------------------------------------------------------------- #
# --- Climatologia mensal e previsão one-step -----------------------------
clim <- read.csv("dados/chuva_climatologia_mensal_sp_1991_2020.csv")
barplot(clim$chuva_media_hist, names.arg = month.name[1:12],
        col = "cornflowerblue", ylab = "mm", main = "Climatologia mensal (1991–2020)")

recente <- read.csv("dados/chuva_mensal_sp_recente.csv")
y <- ts(recente$chuva_mm, start = c(2024, 1), frequency = 12)
plot(y, main = "Chuva mensal recente (2024–2026)")
abline(h = 120, col = "gray50", lty = 2)                   # referência visual

# -------------------------------------------------------------------------- #
# Bloco 26 (extraído da documentação)
# -------------------------------------------------------------------------- #
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
