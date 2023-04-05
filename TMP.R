





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

Allows statements as, for example, $P(RPR> 80\%)$:







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
