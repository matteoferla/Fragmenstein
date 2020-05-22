# Covid and Fragmenstein

Fragmenstein was sparked to life by COVID Moonshot project [the COVID moonshot project](https://discuss.postera.ai/c/covid).
[This dataset](https://github.com/postera-ai/COVID_moonshot_submissions) has some unique peculiarities that potentially
are not encountered in other projects. Namely, humans look at the bound fragments and suggest followups via the form.
However, there are some problems:

* the form does not have a `no inspiration hit` option, so many users submitted `x0072` the first as inspiration when submitting docked libraries.
* the inspiration hits are common to a group of submissions by a user, even if one went into one and another to another.
* some pockets have many hits, so a large amount of not fully overlapping hits are submitted
* some users submitted a mispelt hit

For an example of the script used, see [covid.py](covid.py).

For a comparision of how the three method fair with the daset see [three modes compared](three_modes_compared.md).

## Over-inspired problem

Fragmenstein full-merge mapping works well for two

The 'TRY-UNI-714a760b-1' compound (`Cc1c(N)cncc1NC(=O)CC1CCCCC1`) is purposed to be inspired by x0107, x0434, x0678, x0748, x0995, x1382.

![fragments](images/toomany_frag.png)

This takes forever to make a template... which comes out awful.

![fragments](images/toomany_follow.png)

When placed and minimised the compound drifts off. The reason for this is that there are only two atoms that map.

In reality only x0107, x0678, x0995 were the true inspirations. When this is corrected, the scaffold is basically the followup.

![fragments](images/toomany_perfect.png)

So the question is: how does one fix this?
Before that it is best to see how frequent this is:

![fragments](images/toomany_distro.png)

Of the 2,000 with 1 hit as inspiration, 500 are based upon x0072.
These are not really inspirations, just a case where `null` was not a choice.

The wobbly extras are good to set deviations for the coordinate constraints...