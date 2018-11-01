# Codingame.com

This repository contains the sources - buggy, incomplete and cluttered - in the
authentic form in which I submit them, after the contest is over.

## The Puzzles

The directory containing my submissions to the various _Codingames_ puzzles,
just to showcase how much fun I had. I am actually glad I discovered this, I
find this rather enjoyable to play, more so than video games perhaps.

## The Integer library

A set of reusable code snippets, designed to be copy/pasted in future codingames
puzzles: compact, self-contained, simple to understand and reason about.

### Concept

Since _Codingames_ puzzles play in an arena that is a `x, y` grid of roughly
`-10000,10000`, that allows me to play with fast integers approximate math, such
as approximate square root, approximate sine and cosine calculations, which are
done without branching and converge in a few iterations (hopefully).

### Language

Only presented in C++ for the moment, it contains basic integer vector
operations, with some debbuging functions (we only have print statements
available to help :D) and the aapproximative integer maths.

### Fun With Approximations

Since it's about approximations, within the range of the arena, don't expect
these functions to return good results, in particular when the integers are
small: e.g. computation of the square root of `(1, 1)` will not satisfy the
general use case you are probably looking for.

I quite like this approximate library, because it was built to minimize
overflow/underflow with care and avoid division by 0 using stupid tricks, and so
it forced me to be a little creative, within the confine of the puzzles.

# Licence

The work is not licensed. Anyway it's probably buggy, incomplete, stupid and you
can probably do better in a few hours of work. Don't reuse it. But if you want
to, you can actually.
