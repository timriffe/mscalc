A lightweight, fast, and memory-efficient state occupancy time calculator that supports arbitrary state spaces (as far as I can tell).

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

# and total life expectancy:
sum(Lxs)
```

Also, note, you can provide tidy or wide `data.frames` as `p_list`, in which case they are reformatted internally, but there is an efficiency tradeoff, see:
```
# install.packages("bench")
prep   <- prepare_p_list(transitions_example)
age    <- prep$age
p_list <- prep$p_list

# Tidy version (original)
tidy_input <- transitions_example

# Benchmark
bench::mark(
  list_input = calc_occupancies(p_list = p_list, age = age, origin_state = "P"),
  tidy_input = calc_occupancies(p_list = tidy_input, origin_state = "P"),
  check = TRUE,
  iterations = 100
)
```

That is, if you're really doing big runs (bootstrapping, etc) then it's best to send data in directly in the ideal list format, i.e. a names list of transitions.

Here's a sneak peak at performance in terms of speed:

<img width="560" height="346" alt="image" src="https://github.com/user-attachments/assets/31345f0d-5bd0-424f-9387-630f6ebea429" />

And in terms of memory-usage:

<img width="560" height="346" alt="image" src="https://github.com/user-attachments/assets/54db5e64-430d-4fe8-8d8e-736e7acbc655" />

Check out the performance vignette for more details on these comparisons, and  suggest or prepare others so that we can find out whether this calculator is truly best in some category.

