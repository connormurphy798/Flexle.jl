# Future deprecation notice

Flexle will eventually be supplanted by [WeightVectors](https://github.com/LilithHafner/WeightVectors.jl), a
package using a similar sampling algorithm that solves numerical imprecision problems inherent to Flexle. At
present, `push!` and `deleteat!` are not implemented in WeightVectors, so Flexle may temporarily be preferred
for applications requiring those operations.

# Flexle

Fast, dynamically weighted random sampling.

## Overview

Flexle (<ins>flex</ins>ible, binary-<ins>le</ins>vel rejection sampling) is a package for
**high performance sampling from discrete distributions** under the constraint of being able to
**quickly modify said distributions**.

Other strategies for fast sampling from weighted distributions (such as the
[alias method](https://en.wikipedia.org/wiki/Alias_method)) require precomputing a data structure from which
one can sample quickly but which cannot be "updated"; a change to any of the weights requires recomputing the
entire data structure. Flexle is designed to be dynamic in that, in addition to enabling fast sampling, it
supports fast updating, addition, and removal of weights.

The fundamental sampler data structure and corresponding sampling algorithm are described
[here](https://www.aarondefazio.com/tangentially/?p=58) by Aaron Defazio. Briefly, the strategy involves
grouping weights into "levels" of similar value and tracking the total weight of each level. Sampling is
performed by first using the cumulative distribution function of the levels to select a level, then using
rejection sampling to select an element from the chosen level.

The primary conceptual change to Defazio's method is in allowing **unlimited dynamic range of weights**.
As its weights are updated, a `FlexleSampler` adaptively adds and removes levels to accommodate. This along
with storing weights separately from the levels themselves means that Flexle can in principle be used to model
any discrete distribution, including one where a subset of events happen with probability 0.

For a detailed overview of Flexle's principle of operation (including runtime analysis), see `docs/principle.md`.

## Installation

Flexle is available through Julia's package manager. To install, simply run the following in the Julia REPL:

```
]add Flexle
```

## API

Users interface with Flexle exclusively through instances of the `FlexleSampler` struct. A `FlexleSampler`
is used to store, modify, and sample from a discrete distribution. The `FlexleSampler` is generated from a
collection of `weights` (where `weights[i]` is the weight associated with element `i`) and interacted with using
the following:

| Function | Description |
|---|---|
| `FlexleSampler(weights)` | Create a `FlexleSampler` instance from an `AbstractVector{Float64}` of `weights`. | 
| `getindex(sampler, i)`[^1] | Get the weight of element `i` in `sampler`. |
| `setindex!(sampler, w, i)`[^2] | Update the weight of element  `i` in `sampler` to be equal to `w`, returning the difference between the new and old weights of element `i`. |
| `getweights(sampler)` | Get a `Vector` of all the weights in `sampler`. |
| `numweights(sampler)` | Get the number of weights in `sampler`. |
| `push!(sampler, w)` | Add a new element of weight `w` to `sampler`, placing it at the end of `sampler.weights`. |
| `deleteat!(sampler, i)` | Remove element `i` from `sampler`, shifting the index of every subsequent element over to fill the gap. |
| `sample(sampler)` | Take a single random sample from `sampler`, returning the index of the element sampled. |


[^1]: also `sampler[i]`

[^2]: also `sampler[i] = w`

## Performance

Flexle is a performance-focused package. For a comparison of Flexle's runtime performance to that of
other sampling techniques from StatsBase.jl, see `benchmarking.md`.

TL;DR—In general, Flexle is best suited to use cases where users will _repeatedly_ sample from and modify a
collection of weights, as the linear-time construction of a `FlexleSampler` from such a collection is amortized
by near constant time updating/addition of weights and sampling. That Flexle's sampling is (for most distributions)
asymptotically faster than that of standard linear time algorithms also means that Flexle works best when working
with large collections of weights, where there is a genuine performance improvement to be realized; for smaller
collections, the modest improvement in sampling speed may not be enough to offset the comparatively expensive
sampler initialization.

## Acknowledgements

As described above, I learned of this sampling strategy from Aaron Defazio and his blog post
["Weighted random sampling with replacement with dynamic weights"](https://www.aarondefazio.com/tangentially/?p=58).
Many thanks for the excellent technique and writeup.

This package was originally developed as part of [jOpqua](https://github.com/pablocarderam/jOpqua), an
epidemiological modeling project led by [Pablo Cárdenas R.](https://github.com/pablocarderam) at the Ragon
Institute in Cambridge, MA. Thanks to Pablo for both (1) the motivation for this project, namely the need to
sample and update quickly, and (2) his unending feedback and support thoughout its development.
