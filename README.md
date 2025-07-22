A lightweight and fast state occupancy time calculator that supports arbitrary state spaces (as far as I can tell).

Install using e.g.

```
pak::pak("timriffe/mscalc")
```

or

```
remotes::install_github("timriffe/mscalc")
```

Use like so:

```
library(mscalc)

data(transitions_example)
head(transitions_example)

# impute missing transitions as needed, create list format
p_list <- prepare_p_list(transitions_example)
# calculate occupancy times
occup  <- calc_occupancies(prep$p_list, prep$age, origin_state = "P")
head(occup)

# note:
lxs <- rowSums(occup[,-ncol(occup)])
lx <- rowSums(lxs)
plot(lx)
```
That is, you might want to do some lifetable adjustment to get the Lx equivalent, here we have single-age data:

```
age_interval <- 1
Lxs <- age_interval * (lxs + rbind(lxs[-1,],0)) / 2

# These are total state expectancies
colSums(Lxs)

# and total remaining life expectancy:
sum(Lxs)
```
