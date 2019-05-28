library(tidyverse)
library(rstan)

load("zikachik.Rdata")
tibble(zikachik)


# One island, one disease -----------------

## Select CHIKV in Guadeloupe
GUAD = filter(zikachik,ISLAND=="GUADELOUPE",VIRUS=="CHIKV",WEEK>=5)
## Format data for rstan
GUAD_L = list(W=dim(GUAD)[[1]],
              O_t=GUAD$NCASES,
              Ostar_t=GUAD$Ostar,
              sumO_t=GUAD$CUM_NCASES,
              pop=GUAD$POP[[1]])
## Sample
S_GUAD = stan("TSIR_one_island.stan",data=GUAD_L)
## Results
print(S_GUAD,pars=c("beta","rho","phi"))
## Plot
pred = summary(S_GUAD,pars="lp")[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  mutate(WEEK=GUAD$WEEK)
ggplot() +
  geom_ribbon(data=pred,aes(x=WEEK,ymin=X2.5.,ymax=X97.5.),alpha=.5,fill="grey70") +
  geom_line(data=pred,aes(x=WEEK,y=X50.)) +
  geom_point(data=GUAD,aes(x=WEEK,y=NCASES),shape=21,fill="tomato") +
  labs(x="Week",y="N") +
  theme_bw()
ggsave(file="Figures/fit.png",width=6,height=4)
