Rufibach et al. (2016) give a closed a formula for the distribution of $RPR:=P_{\Delta}(Z\leq - z_{1-\alpha})$, where $\Delta \sim N(d,\tilde\sigma^2/n_0)$, and discuss the shape under different prior choices (RPR='Random probability to reject`, see Kunzmann et al.)

```{r}
#| echo: false

x <- seq(0.001,0.999,0.001)

n <- 100
delta_star <- 0.035
prior_mean <- 0
prior_mean_skeptical <- delta_star
alpha <- 0.05

n_0 <- 6.6

## Rufibach 2016: formula (4)
y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- data.frame(x, y, type="Enthusiastic")

n_0 <- 6.6

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_skeptical-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="Skeptical"))

n_0 <- 25

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="Informative"))

n_0 <- 0.5

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="Noninformative"))

data_output2 <- data_output %>% select(type, x=AP) %>% mutate(y=0)

ggplot(data_power, aes(x, y, colour=type))+geom_line(linewidth=1)+geom_point(data=data_output2, aes(x=x, y=y, colour=type))+scale_colour_manual("Prior type", values=c("purple", "orange", "blue", "darkred"))+theme_bw()+theme(panel.grid = element_blank())+ylab("Density")+xlab("RPR")+ylim(c(0,2))
```





### Hybrid approach: Prior on $p_1$ and $p_0$

A hybrid approach assumes that the true treatment effect $\delta$ is a realization from a random variable $\Delta$ with prior probability density function $p(\delta)$. The prior is defined before data collection ('design prior'). If we assume that $\pi_i:=p_i(\delta)=Beta(a_i,b_i)$, then the posterior distribution is $Beta(a_i+\sum_{k\leq n_i}Y_{i,k},b_i+n_i-\sum_{k\leq n_i}Y_{i,k})$. The variance is given as

$$
  \frac{\left(a_i+\sum_{k\leq n_i}Y_{i,k}\right)\left(b_i+n_i-\sum_{k\leq n_i}Y_{i,k}\right)}{(a_i+b_i+n_i)^2(a_i+b_i+n_i+1)}
$$
  
  ```{r}
#| echo: false
#| results: hide
#| 
a_0 <- 0.01
b_0 <- (1-p_0)/p_0*a_0
a_0/(a_0+b_0)
a_0+b_0
a_1 <- 0.01
b_1 <- (1-p_1)/p_1*a_1
a_1/(a_1+b_1)

## Uniform
a_0 <- 1
b_0 <- 1
a_0/(a_0+b_0)
a_0+b_0
a_1 <- 1
b_1 <- 1
a_1/(a_1+b_1)


## Uniform
a_0 <- 100
b_0 <- 1900
a_0/(a_0+b_0)
a_0+b_0
a_1 <- 100
b_1 <- 1900
a_1/(a_1+b_1)

a_1+b_1

x <- seq(0,1,0.01)
prior_0 <- dbeta(x, a_0, b_0)
prior_1 <- dbeta(x, a_1, b_1)

posterior_0 <- dbeta(x, a_0+n_0*p_0, b_0+n_0-n_0*p_0)
posterior_1 <- dbeta(x, a_1+n_1*p_1, b_1+n_1-n_1*p_1)

var_posterior_0 <- ((a_0+n_0*p_0)*(b_0+n_0-n_0*p_0))/((a_0+b_0+n_0)^2*(a_0+b_0+n_0+1))
var_posterior_1 <- ((a_1+n_1*p_1)*(b_1+n_1-n_1*p_1))/((a_1+b_1+n_1)^2*(a_1+b_1+n_1+1))
sqrt(var_posterior_0+var_posterior_1)


delta_range <- seq(from=-0.1, to=0.1, 0.001)

data_output <- c()

prior_mean <- 0

data_output <- rbind(data_output, data.frame(delta_range, type="Enthusiastic", y=dnorm(delta_range, mean=prior_mean, sd=sqrt(var_posterior_0+var_posterior_1))))

prior_mean <- delta_star

data_output <- rbind(data_output, data.frame(delta_range, type="Skeptical", y=dnorm(delta_range, mean=prior_mean, sd=sqrt(var_posterior_0+var_posterior_1))))

ggplot(data_output, aes(x=delta_range, y=y, colour=factor(type)))+geom_line(linewidth=1.2)+theme_bw()+xlab("delta=p1-p0")+ylab("Density")+scale_x_continuous(breaks=c(-0.1, 0, 0.035, 0.1))+theme(panel.grid = element_blank())+geom_vline(xintercept=delta_star, linetype="dotted")+geom_vline(xintercept=0, linetype="dashed")+scale_linetype("Prior type") + stat_function(fun = dnorm, args = list(mean = prior_mean_enthusiastic, sd = sd_tilde/n_00), xlim = c(0.035,0.1), geom = "area", show.legend = F, fill="purple", alpha=0.3, colour=NA) + stat_function(fun = dnorm, args = list(mean = prior_mean_skeptical, sd = sd_tilde/n_00), xlim = c(-0.1,0), geom = "area", show.legend = F, fill="darkred", alpha=0.3, colour=NA) +scale_colour_manual("Prior type", values=c("purple", "orange", "blue", "darkred"))
```



## Implications


```{r}
#| echo: false

x <- seq(0.001,0.999,0.001)

n <- 100
delta_star <- 0.035
prior_mean_enthusiastic <- 0
prior_mean_skeptical <- delta_star
alpha <- 0.05

n_0 <- 6.6

## Rufibach 2016: formula (4)
y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- data.frame(x, y, type="Enthusiastic")

n_0 <- 6.6

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_skeptical-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="Skeptical"))

n_0 <- 25

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="Informative"))

n_0 <- 0.5

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="Noninformative"))

data_output2 <- data_output %>% select(type, x=AP) %>% mutate(y=0)

ggplot(data_power, aes(x, y, colour=type))+geom_line(linewidth=1)+geom_point(data=data_output2, aes(x=x, y=y, colour=type))+scale_colour_manual("Prior type", values=c("purple", "orange", "blue", "darkred"))+theme_bw()+theme(panel.grid = element_blank())+ylab("Density")+xlab("RPR")+ylim(c(0,2))
```

## Implications

- unimodal if $n=n_0$
- 'hill'-shape if $n_0>n$ (rather unrealistic interpretation)

```{r}
#| echo: false

x <- seq(0.001,0.999,0.001)

n_0 <- 100
n <- 100
delta_star <- 0.035
prior_mean_enthusiastic <- 0
alpha <- 0.05

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- data.frame(x, y, type="n0=100 (enthusiastic)")

n_0 <- 1000

y <- sqrt(n_0/n)*dnorm(-sqrt(n_0/n)*qnorm(1-alpha)-sqrt(n_0)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(n_0/n)*qnorm(x))*(dnorm(qnorm(x)))^(-1)

data_power <- rbind(data_power, data.frame(x, y, type="n0=1000 (enthusiastic)"))

ggplot(data_power, aes(x, y, colour=type))+geom_line(linewidth=1)+scale_colour_manual("Prior sample size", values=c("blue", "turquoise"))+theme_bw()+theme(panel.grid = element_blank())+ylab("Density")+xlab("RPR")+ylim(c(0,5))
```

## Implications

Allows statements as, for example, $P(RPR> 80\%)$:

```{r}
#| echo: false

x <- 0.8

prior_mean_enthusiastic <- 0

n_0_range <- c(0.5, 6.6, 25, 50, 100)
n_range <- seq(0.01, 1000, by=1)

data_rpr <- expand.grid(n_0_range, n_range)
names(data_rpr) <- c("n_0_range", "n_range")
data_rpr$p_rpr <- pnorm(-sqrt(data_rpr$n_0_range/data_rpr$n_range)*qnorm(1-alpha)-sqrt(data_rpr$n_0_range)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(data_rpr$n_0_range/data_rpr$n_range))
data_rpr$type <- "Enthusiastic"

prior_mean_enthusiastic <- delta_star

n_0_range <- c(0.5, 6.6, 25, 50, 100)
n_range <- seq(0.01, 1000, by=1)

data_rpr_skep <- expand.grid(n_0_range, n_range)
names(data_rpr_skep) <- c("n_0_range", "n_range")
data_rpr_skep$p_rpr <- pnorm(-sqrt(data_rpr_skep$n_0_range/data_rpr_skep$n_range)*qnorm(1-alpha)-sqrt(data_rpr_skep$n_0_range)/sd_tilde*(prior_mean_enthusiastic-delta_star)-sqrt(data_rpr_skep$n_0_range/data_rpr_skep$n_range))
data_rpr_skep$type <- "Skeptical"

data_rpr <- rbind(data_rpr, data_rpr_skep)

ggplot(data_rpr, aes(x=n_range, y=p_rpr, colour=factor(n_0_range), linetype=type))+geom_line(linewidth=1)+theme_bw()+theme(panel.grid = element_blank())+scale_color_brewer("Prior sample size", palette="Reds")+xlab("Final sample size (one arm)")+geom_hline(yintercept=0.8, linetype="dashed")+ylab("P(RPR>0.8)")+scale_x_continuous(breaks=c(0, 100, 250, round(min(data_rpr$n_range[data_rpr$n_0_range==50 & data_rpr$p_rpr>=0.8]),0), 500, 750, 1000))+geom_vline(xintercept=min(data_rpr$n_range[data_rpr$n_0_range==50 & data_rpr$p_rpr>=0.8]), linetype="dashed")+geom_vline(xintercept=100, linetype="dashed")+scale_y_continuous(breaks=c(0,0.5,0.8,1))+scale_linetype("Prior type")
```

## Decomposition of AP

- AP integrates over whole $\delta$ range, including 'non-favorable'
$$
AP=\overbrace{P(Z\leq - z_{1-\alpha}, \Delta>\delta^*)}^{(1)}+\overbrace{P(Z\leq - z_{1-\alpha}, 0<\Delta\leq\delta^*)}^{(2)} \\
+\overbrace{P(Z\leq - z_{1-\alpha}, \Delta\leq0)}^{(3)}
$$
- AP=(1) Probability of Type-I error + (2) 'Non-inferior, but treatment effect not relevant' + (3) 'Non-inferior, treatment effect relevant'
- $AP \approx P(Z\leq - z_{1-\alpha}, \Delta\leq0)$ (see Spiegelhalter), pharma (favors AP, shortterm) vs regulators (favors the other, longterm), see discussion in Kunzmann et al.

## Decomposition of AP

Under enthusiastic prior:

```{r}
#| echo: false

library(mvtnorm)
library(tidyverse)
library(ggplot2)

delta_star <- 0.035
prior_mean_enthusiastic <- 0
n_0 <- 6.6
n <- 100
sigma_sim <- matrix(c((n_0+n)/n, 1, 1, 1), nrow=2, ncol=2, byrow=T)

data_sim <- data.frame(rmvnorm(n=100000, mean=c(prior_mean_enthusiastic, prior_mean_enthusiastic), sd_tilde^2/n_0*sigma_sim))
names(data_sim) <- c("delta_pred", "delta")

data_sim$region <- ifelse(data_sim$delta<=0 & data_sim$delta_pred<=-qnorm(1-alpha)*sd_tilde/sqrt(n)-(prior_mean_enthusiastic-delta_star), 1, 0)
data_sim$region <- ifelse(data_sim$delta>0 & data_sim$delta<=delta_star & data_sim$delta_pred<=-qnorm(1-alpha)*sd_tilde/sqrt(n)-(prior_mean_enthusiastic-delta_star), 2, data_sim$region)
data_sim$region <- ifelse(data_sim$delta>delta_star & data_sim$delta_pred<=-qnorm(1-alpha)*sd_tilde/sqrt(n)-(prior_mean_enthusiastic-delta_star), 3, data_sim$region)

data_sim$region <- factor(data_sim$region, levels=0:3, labels=c("'Not significant'", "(3)", "(2)", "(1)"))

ggplot(data_sim, aes(x=delta, y=delta_pred, colour=factor(region)))+geom_point(size=0.3)+theme_bw()+theme(panel.grid.minor = element_blank())+geom_vline(xintercept=0.035, linetype="dashed")+geom_hline(yintercept=-qnorm(1-alpha)*sd_tilde/sqrt(n)-(prior_mean_enthusiastic-delta_star), linetype="dotted")+labs(caption="Dotted line: Upper confidence interval from mean of predicted delta <= non-inferiority margin")+ylab("Predictive distribution of Delta")+xlab("Delta")+scale_color_brewer("Region", palette="Reds")

data.frame(data_sim %>% group_by(region) %>% summarise(prop=n()/nrow(data_sim)))
```

## Expected power

$$
AP=\overbrace{P(Z\leq - z_{1-\alpha}, \Delta>\delta^*)}^{(1)}+\overbrace{P(Z\leq - z_{1-\alpha}, 0<\Delta\leq\delta^*)}^{(2)} \\
+\overbrace{P(Z\leq - z_{1-\alpha}, \Delta\leq0)}^{(3)}
$$
For our noninferiority setting consider (2)+(3).

$$
P(Z\leq - z_{1-\alpha}, \Delta\leq\delta^*)=P(Z\leq - z_{1-\alpha}|\Delta\leq\delta^*)P(\Delta\leq\delta^*)\\=\underbrace{E\left[P_{\Delta\leq\delta^*}(Z\leq - z_{1-\alpha})\right]}_{EP}P(\Delta\leq\delta^*),
$$
where $EP$ is **'expected power'** (Kunzmann et al.).

## To make things more powerful confusing

Spiegelhalter calls $P(Z\leq - z_{1-\alpha}, \Delta\leq\delta^*)$ the **'prior adjusted power'** (PAP):

$$
\underbrace{P(Z\leq - z_{1-\alpha}, \Delta\leq\delta^*)}_{PAP}=\underbrace{E\left[P_{\Delta\leq\delta^*}(Z\leq - z_{1-\alpha})\right]}_{EP}\underbrace{P(\Delta\leq\delta^*)}_{constant}
$$

## EP, PAP, AP

Noninferiority setting $\Delta\leq\delta^*$:

```{r}
#| echo: false

n <- 100
delta_star <- 0.035
prior_mean_skeptical <- delta_star
prior_mean_enthusiastic <- 0


set.seed(1)
draws <- rnorm(100000, mean=prior_mean_skeptical, sd=sd_tilde/sqrt(6.6))

power_classic <- pnorm(-qnorm(1-alpha)-sqrt(n)/sd_tilde*(draws-delta_star))

data_ep <- data.frame(ep=mean(power_classic[draws<=delta_star]), pap=mean(power_classic[draws<=delta_star])*pnorm(delta_star, mean=prior_mean_skeptical, sd=sd_tilde/sqrt(6.6)), ap=mean(power_classic), const=pnorm(delta_star, mean=prior_mean_skeptical, sd=sd_tilde/sqrt(6.6)), type="Skeptical")

set.seed(1)
draws <- rnorm(100000, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(6.6))

power_classic <- pnorm(-qnorm(1-alpha)-sqrt(n)/sd_tilde*(draws-delta_star))

data_ep <- rbind(data_ep, data.frame(ep=mean(power_classic[draws<=delta_star]), pap=mean(power_classic[draws<=delta_star])*pnorm(delta_star, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(6.6)), ap=mean(power_classic), const=pnorm(delta_star, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(6.6)), type="Enthusiastic"))

set.seed(1)
draws <- rnorm(100000, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(25))

power_classic <- pnorm(-qnorm(1-alpha)-sqrt(n)/sd_tilde*(draws-delta_star))

data_ep <- rbind(data_ep, data.frame(ep=mean(power_classic[draws<=delta_star]), pap=mean(power_classic[draws<=delta_star])*pnorm(delta_star, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(25)), ap=mean(power_classic), const=pnorm(delta_star, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(25)), type="Informative"))

set.seed(1)
draws <- rnorm(100000, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(0.5))

power_classic <- pnorm(-qnorm(1-alpha)-sqrt(n)/sd_tilde*(draws-delta_star))

data_ep <- rbind(data_ep, data.frame(ep=mean(power_classic[draws<=delta_star]), pap=mean(power_classic[draws<=delta_star])*pnorm(delta_star, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(0.5)), ap=mean(power_classic), const=pnorm(delta_star, mean=prior_mean_enthusiastic, sd=sd_tilde/sqrt(0.5)), type="Noninformative"))


data_ep$ep <- round(data_ep$ep, 4)
data_ep$ap <- round(data_ep$ap, 4)
data_ep$pap <- round(data_ep$pap, 4)
data_ep$const <- round(data_ep$const, 4)

data_ep
```

## Conditional Bayesian power

Uses posterior distribution to define 'Bayesian significance': $S^B:=P(\Delta\leq\delta^*|data)=1-\epsilon$, (not $\alpha$)
- Remember $D=\overline{Y}_{1}-\overline{Y}_{0}$ with $D \sim N(\delta, \frac{\tilde\sigma^2}{n})$, where $\tilde\sigma^2=\sigma_1^2+\sigma_1^2$
Suppose that $\Delta \sim N(d, \tilde\sigma^2/n_0)$ then the posterior distribution given $D$ is distributed as $\Delta|D \sim N\left(\frac{n_0d+nD}{n_0+n}, \frac{\tilde\sigma^2}{n_0+n}\right)$. We are interested in

$$

$$

- $P(S^B|\delta)=\Phi\left(-\sqrt{\frac{n_0+n}{n}}z_{1-\epsilon}-\frac{\sqrt{n}}{\tilde\sigma}\left[\delta-\frac{n_0+n}{n}\delta^*+\frac{n_0}{n}d\right]\right)$
- $P(S^{F}|\delta)=\Phi\left(-z_{1-\alpha}-\frac{\sqrt{n}}{\tilde\sigma}(\delta-\delta^*)\right)$

## Conditional Bayesian power

```{r}
#| echo: false

n <- 100
delta_star <- 0.035


prior_mean <- 0
n_0 <- 6.6

delta <- seq(-0.025, 0.05, 0.001)

y_bayes <- pnorm(-(delta)*sqrt(n)/sd_tilde+(delta_star)*(n_0+n)/(sd_tilde*sqrt(n))-(prior_mean*n_0)/(sd_tilde*sqrt(n))-sqrt((n_0+n)/n)*qnorm(1-alpha))

data_power <- data.frame(delta, y=y_bayes, type="Bayes: Enthusiastic")

prior_mean <- 0.035
n_0 <- 6.6

y_bayes <- pnorm(-(delta)*sqrt(n)/sd_tilde+(delta_star)*(n_0+n)/(sd_tilde*sqrt(n))-(prior_mean*n_0)/(sd_tilde*sqrt(n))-sqrt((n_0+n)/n)*qnorm(1-alpha))
data_power <- rbind(data_power, data.frame(delta, y=y_bayes, type="Bayes: Skeptical"))


y_freq <- pnorm(-(delta-delta_star)*sqrt(n)/sd_tilde-qnorm(1-alpha))
data_power <- rbind(data_power, data.frame(delta, y=y_freq, type="Frequentist"))

ggplot(data_power, aes(x=delta, y=y, colour=type))+geom_line()+theme_bw()+theme(panel.grid.minor = element_blank())+scale_colour_brewer("", palette="Dark2")+ylab("Probability of rejection")+ylab("Probability to reject")+scale_y_continuous(breaks=c(0,0.5,0.8,1))+ggtitle("n=100 (one arm)")+geom_vline(xintercept=0, linetype="dashed")+geom_vline(xintercept=0.035, linetype="dashed")

data_power$y <- round(data_power$y, 2)
data_power %>% filter(delta%in%c("0", "0.035"))
```

## Average Bayesian power

- $P(S^B)=\Phi\left(-\sqrt{\frac{n_0}{n}}z_{1-\epsilon}-\frac{\sqrt{n_0}\sqrt{n_0+n}(d-\delta^*)}{\sqrt{n}\tilde\sigma }\right)$

```{r}
#| echo: false

### Hybrid
n <- 100

n_0 <- c(6.6)

data_output <- data.frame(type="Enthusiastic", n_0, n, AP=round(pnorm(sqrt(n_0/(n_0+n))*(-qnorm(1-alpha)-sqrt(n/(sd_tilde^2))*(prior_mean_enthusiastic-delta_star))),3), upper_AP=round(pnorm(-(prior_mean_enthusiastic-delta_star)*sqrt(n_0/sd_tilde^2)), 3))

n_0 <- c(6.6)

data_output <- rbind(data_output, data.frame(type="Skeptical", n_0, n, AP=round(pnorm(sqrt(n_0/(n_0+n))*(-qnorm(1-alpha)-sqrt(n/(sd_tilde^2))*(prior_mean_skeptical-delta_star))),3), upper_AP=round(pnorm(-(prior_mean_skeptical-delta_star)*sqrt(n_0/sd_tilde^2)), 3)))

n_0 <- c(25)

data_output <- rbind(data_output, data.frame(type="Informative", n_0, n, AP=round(pnorm(sqrt(n_0/(n_0+n))*(-qnorm(1-alpha)-sqrt(n/(sd_tilde^2))*(prior_mean_enthusiastic-delta_star))),3), upper_AP=round(pnorm(-(prior_mean_enthusiastic-delta_star)*sqrt(n_0/sd_tilde^2)), 3)))

n_0 <- c(0.5)

data_output <- rbind(data_output, data.frame(type="Noninformative", n_0, n, AP=round(pnorm(sqrt(n_0/(n_0+n))*(-qnorm(1-alpha)-sqrt(n/(sd_tilde^2))*(prior_mean_enthusiastic-delta_star))),3), upper_AP=round(pnorm(-(prior_mean_enthusiastic-delta_star)*sqrt(n_0/sd_tilde^2)), 3)))

data_output <- data_output %>% select(type, AP)

n <- 100
delta_star <- 0.035


prior_mean <- 0
n_0 <- 6.6

y <- pnorm(-sqrt(n_0/n)*qnorm(1-alpha)-((prior_mean-delta_star)*sqrt(n_0+n)*sqrt(n_0))/(sd_tilde*sqrt(n)))

data_output2 <- data.frame(type="Enthusiastic", n_0, n, AP=round(y, 3))

prior_mean <- 0.035
n_0 <- 6.6

y <- pnorm(-sqrt(n_0/n)*qnorm(1-alpha)-((prior_mean-delta_star)*sqrt(n_0+n)*sqrt(n_0))/(sd_tilde*sqrt(n)))

data_output2 <- rbind(data_output2, data.frame(type="Skeptical", n_0, n, AP=round(y, 3)))

prior_mean <- 0
n_0 <- 25

y <- pnorm(-sqrt(n_0/n)*qnorm(1-alpha)-((prior_mean-delta_star)*sqrt(n_0+n)*sqrt(n_0))/(sd_tilde*sqrt(n)))

data_output2 <- rbind(data_output2, data.frame(type="Informative", n_0, n, AP=round(y, 3)))

prior_mean <- 0
n_0 <- 0.5

y <- pnorm(-sqrt(n_0/n)*qnorm(1-alpha)-((prior_mean-delta_star)*sqrt(n_0+n)*sqrt(n_0))/(sd_tilde*sqrt(n)))

data_output2 <- rbind(data_output2, data.frame(type="Noninformative", n_0, n, AP=round(y, 3)))

data_output2$AP_bayes <- data_output2$AP

data_output <- left_join(data_output, data_output2 %>% select(type, AP_bayes), by="type")

data_output
```

Of course, Bayesian average power can again be decomposed into EP and PAP.

## Power vocabulary summary (`r emoji("biceps")`)

- 'Power': Classical power which is **conditional** on a fixed alternative point estimate. Probability to reject.
- Bayesian power: Using posterior distribution.
- Conditional power: Frequentist power in interim analysis (not discussed here).
- Average power: **Marginal** probability of rejection, classical power 'averaged' over (design) prior. Assurance. Bayesian **predictive** power.
- Expected power: Weighted average of the probability to reject in a 'relevant' region. Also called conditional expected power (`r emoji("confused")`).
- Probability of success (PoS): Very often authors use PoS for 'average power'. Depends on success definition.

## From the SAFE-SSPE study protocol

For sample size calulcation:

*We used a Monte-Carlo simulation approach based on an Agresti-Caffo confidence interval for risk difference because normal approximation formulas do not hold in case of rare events*