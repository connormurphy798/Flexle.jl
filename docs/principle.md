# Principle of operation

To enable fast random sampling from a modifiable collection of weights, Flexle implements
the data structure and sampling algorithm described by Aaron Defazio in his 2016 blog post,
["Weighted random sampling with replacement with dynamic weights"](https://www.aarondefazio.com/tangentially/?p=58).[^terminology]
What follows is a discussion of how Flexle works and the runtime of its methods.

Given a collection of weights $W$ where $W_i$ corresponds to the weight associated with element
i, Flexle creates a `FlexleSampler` object that organizes said weights into "levels" according
to magnitude. Sampling is performed by first selecting a level according to the cumulative
distribution function of the levels (CDF sampling), then selecting an element from within the
chosen level using rejection sampling.

## Sampling data structure

A `FlexleSampler` is made up of levels, each of which is contains all elements in the
sampler with weights between two adjacent powers of 2, e.g. in $[8, 16)$. This choice of
bounds has two benefits: [^fastlog2]

1. The maximum weight in a level is strictly less than twice the minimum weight, which makes
rejection sampling fast.
2. Levels being logarithmic in scale allows for efficient representation of wide dynamic ranges
of weights.

The number and range of levels is dynamically chosen and adjusted based on the minimum and maximum
nonzero weights in a `FlexleSampler`. (Specifically, a sampler made from a vector of weights $W$
will have $\lceil \log_2(\frac{\max(W^+)}{\min(W^+)}) \rceil + 1$ levels, where $W^+$ is the set
of nonzero weights in $W$.) For example, `s = FlexleSampler([2.0, 1.5, 2.5, 0.0, 0.3, 3.5])`
is represented as follows:

| Level | Range | Elements |
| - | - | - |
| $1$ | $[2.0, 4.0)$ | `[1, 3, 6]` |
| $2$ | $[1.0, 2.0)$ | `[2]` |
| $3$ | $[0.5, 1.0)$ | `[]` |
| $4$ | $[0.25, 0.5)$ | `[5]` |

Each level also tracks both the sum of its weights and its maximum weight for use in the sampling
algorithm (details [below](#step-1:-cdf-sample-a-level)).

| Level | Range | Elements | Sum | Max |
| - | - | - | - | - |
| $1$ | $[2.0, 4.0)$ | `[1, 3, 6]` | `8.0` | `3.5` |
| $2$ | $[1.0, 2.0)$ | `[2]` | `1.5` | `1.5` |
| $3$ | $[0.5, 1.0)$ | `[]` | `0.0` | `0.0`[^max] | 
| $4$ | $[0.25, 0.5)$ | `[5]` | `0.3` | `0.3` |

Finally, the sampler tracks the sum of all its weights (or equivalently, the sum of the sums of
its levels); in this case, `s.sum = 9.8`.

Given a vector $W$ of $n$ weights, initializing a `FlexleSampler` takes $Θ(n + \log(\frac{\max(W^+)}{\min(W^+)})$
time; the second term—the $\log$ of the dynamic range of the weights—gives the time to create the
levels, and the first gives the time to populate them. While the the number of weights and the
$\log$ of their dynamic range are not formally related, for common distributions the former tends
to rapidly outpace the later as $n$ grows large, making `FlexleSampler` initialization generally
linear in the number of weights.

## Sampling algorithm

The user can draw a single random sample from a `FlexleSampler` using `sample(::FlexleSampler)`;
the function returns the index of the element sampled. The algorithm for sampling is in two steps:
CDF sampling a level, then rejection sampling within said level.

### Step 1: CDF sample a level

For a general ordered collection of $n$ weights $V$ with sum $s$, CDF sampling is a sampling
algorithm whose runtime is linear in the number of weights. Its operation is as follows:

1. Select a random number $r \in [0, s)$.
2. Initialize a cumulative sum $c$ at $0$.
3. For $i \in [1 .. n]$, add $V_i$ to $c$.
4. Return the first $i$ such that $c \geq r$.

`Flexle.sample`'s first step is to select a level from a `FlexleSampler` using CDF sampling.
The "weight" corresponding to each level is equal to the sum of the weights of its elements;
this value is precomputed as part of the data structure. Because a `FlexleSampler` with weights
$W$ contains $\lceil \log_2(\frac{\max(W^+)}{\min(W^+)}) \rceil + 1$ levels, this portion of the
algorithm runs in $O(\log(\frac{\max(W^+)}{\min(W^+)}))$ time; in other words, time to sample is
proportional to the $\log$ of the dynamic range of the weights.

### Step 2: Rejection sample an element from the level

For a general ordered collection of $n$ weights $V$ with maximum weight $m$, rejection sampling operates
as follows:

1. Select a random $i \in [1 .. n]$.
2. Select a random $r \in [0, m)$.
3. If $V_i >= r$, return $i$. Else, return to step 1.

Once `Flexle.sample` has selected a level, the next step is to use rejection sampling on the
elements in said level, taking advantage of the precomputed max of any weight in the level;
the element selected through rejection sampling is that returned by `Flexle.sample`. 

The runtime of rejection sampling depends on the number of iterations $T$, which is itself proportional
to the average probability of "rejecting" a choice of sample.

$$E(T) = \sum_{t=1}^{\infty}{t \cdot P_\text{avg}(\text{reject})}^t$$

The probability of rejection for a choice of
sample is $P(\text{reject}_i) = 1 - \frac{V_i}{m}$. The average probability of rejection is then:

$$P_\text{avg}(\text{reject}) = \frac{1}{n}\sum_{i=1}^{n}{(1 - \frac{V_i}{m})}$$

Rejection sampling is therefore most efficient when the above value is minimized. By construction,
the maximum weight in a `FlexleSampler` level is strictly less than double the minimum. The worst
possible case for within-level rejection sampling, then, is when all of the following hold:

- the number of elements in the level $n$ approaches infinity
- one element (arbitrarily, say element 1) has a weight which approaches the level's upper bound $u$
- every other element has a weight equal to the lower bound $l$

In such a case, the average probability of rejection is given by:

$$P_\text{avg}(\text{reject}) = \lim_{n \to \infty}[\frac{1}{n}\sum_{i=1}^{n}{(1 - \frac{V_i}{m})}] = (1 - \frac{V_1}{m}) + \lim_{n \to \infty}[\frac{1}{n}\sum_{i=2}^{n}{(1 - \frac{V_i}{m})}] = 1 - \frac{m}{m} + \frac{l}{u} = \frac{1}{2}$$

The expected value of iterations in this case is therefore:

$$E(T) = \sum_{t=1}^{\infty}{(t \cdot (\frac{1}{2}})^t) = 2$$

In other words, the expected number of rejection sampling iterations has a constant upper bound.
Within-level rejection sampling in a `FlexleSampler` therefore runs in expected $\Theta(1)$ time.

### Summary



[^terminology]: Defazio does not assign names to the data structure or algorithm he describes. For
our purposes, we have named the data structure a "flexle sampler" (short for "flexible, binary-level
sampler; implemented by the `FlexleSampler` struct in the Flexle package) and the corresponding
algorithm "flexle sampling" (implemented by the `sample(::FlexleSampler)` method).

[^fastlog2]: A third, minor benefit is that using powers of $2$ as bounds enables fast calculation
of the level in which a weight belongs using its IEEE 754 floating point representation; however, this
is not essential to implementing the flexle sampler data structure efficiently.

[^max]: Strictly speaking, the maximum weight is undefined for levels containing no weights.
