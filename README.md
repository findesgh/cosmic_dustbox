# A toolbox for cosmic dust.
In 1930, after confirming that distant stars in the Milky Way are dimmed more
strongly than one would expect from geometry, Trumpler conjectured the
existence of "dust particles of various sizes" in the interstellar space of our
Galaxy. Countless observations have since then confirmed this and the
"interstellar extinction" remains one of the observational pillars of the study
of cosmic dust. It was later discovered that there is dust also in
interplanetary space and that it is important for planet formation, there is
dust in other galaxies, even far away galaxies, and maybe there also is dust in
intergalactic space. It was confirmed that dust grains do indeed come in
different sizes and, moreover, in different shapes and different chemical
compositions. In short, while many things have been understood, it also became
apparent that cosmic dust is very complicated. The goal of this software
package is to make work with models of cosmic dust on the computer easier while
still providing control over all the involved intricacies. An ambitious goal
for sure, but we'll try our best.

To introduce the objects central to how we describe cosmic dust, assume we are
interested in the dust in a certain environment. Assume then that we cut some
volume out of said environment. The volume should, of course, be representative
the environment. So if we are interested in the dust in the Solar System as a
whole, a cubic meter of space somewhere between Earth and Mars will be too
small and a box the size of the Milky Way will be too large. Suppose we
collected all of the dust grains contained in this volume. We would like to be
able to describe this collection with a reasonable degree of accuracy without
having to store information on every single grain. If we investigated the
properties of single grains, at first the diversity would probably be
overwhelming since hardly any two grains would be exactly alike. But after
looking at enough grains, we would start discovering similarities. Grains that
have similar structure, grains that look like larger or smaller versions of
other grains, grains that have similar chemical composition and so on. So if we
are fine with some level of abstraction and generalization, our grain
collection could be represented like this:

![](figs/grain-pop.png)

with color-coded chemical composition. How do we go about categorizing these
grains so that we can think in these categories without worrying about
individual grains? Two obvious properties to group by are the chemical
composition and the shape:

![](figs/grain-species.png)

The "blueprint grain", which we can rescale to obtain any of the grains in one
of the groups, is referred to as a grain species. So one of our groups is
simply a grain species combined with a size distribution and we call it a grain
population. If we think back to the abstract representation of our collection
of grains, it is not very different from one of the groups. It simply has more
than one grain species, each with its associated size distribution. So to keep
it simple, we will also call this a grain population.

So, to recap, what we need to describe a collection of grains are:

- grain species
- size distributions
- grain populations

This approach can be refined as needed by simply splitting grain species. In
the extreme, we would have one species per grain and each species would have a
delta-peak size distribution and we would be back to the overwhelming diversity
we started out from.

# Literature
"Physics of the Interstellar and Intergalactic Medium" by B. Draine provides
a good introduction (and more than that) to cosmic dust.
