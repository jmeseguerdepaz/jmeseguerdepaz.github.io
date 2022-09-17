---
layout: post
title:  "Signed distance function for NACA00xx"
date:   2021-04-12 19:30:00 +0100
tags: [graphics, math, sdf, shadertoy]
math: true
img_path: /img/2021-04-12-signed distance function for NACA00x
description: 'For a recent hobby fluid simulator project I wanted to use a NACA airfoil. This post explains how I went about creating it.'
---

<iframe width="640" height="360" frameborder="0" src="https://www.shadertoy.com/embed/ssBGWd?gui=true&t=10&paused=false&muted=false" allowfullscreen></iframe>

# Introduction

For a recent hobby fluid simulator project I wanted to use a NACA airfoil. The solver already handled analytic signed distance fields to determine boundaries, so all I needed was one for the NACA. However, I didn't find anyone readily available. This post explains how I went about creating it. This was the final result:

# Signed distance function
A [signed distance function](https://en.wikipedia.org/wiki/Signed_distance_function "signed distance function") is a function that returns the signed distance from a given input point (let's call it $p=[p_x, p_y]^T$) to a given target (in our case, a NACA00xx). The signed distance is simply the distance from the point to the target as a positive value if the input point is outside the target, or a as negative value if the input point is within the target.

This leads us to the concept of *signed distance field*, which simply the field resulting from evaluating a signed distance function at every point in space.

Signed distance functions have many applications, because they are a compact way of representing geometries that would be very difficult to represent otherwise (such as fractals, fluids of complex topology, etc). This make them very attractive for applications such as computer graphics and simulation.

The best references for signed distance functions are, most likely, Iñigo Quilez's articles on them:
- [distance functions](https://iquilezles.untergrund.net/www/articles/distfunctions/distfunctions.htm)
- [2D distance functions](https://iquilezles.untergrund.net/www/articles/distfunctions2d/distfunctions2d.htm)

# NACA00xx
The [NACA00xx](https://en.wikipedia.org/wiki/NACA_airfoil "NACA00xx") is a family of airfoil shapes as defined by the National Advisory Committee for Aeronautics (NACA).

There are many series of NACA airfoils, some defined by two digits (NACA00xx), other defined by four, five or more.

The 00xx series is a family of symmetric airfoils that are defined by their thickness (expressed in the last two digits, the xx, as a percentage of the chord -the length of the airfoil-).

The airfoil shape is given by a function:
\\[
y_t=5.0 t \left(0.2969 \sqrt{x} - 0.1037 x^{4} + 0.2843 x^{3} - 0.3516 x^{2} - 0.126 x\right)
\\]
where:
$x \in [x,1]$ is the chord coordinate.
$t \in [0, 0.4]$ is the thickness as a percentage of the chord length.

To get the full profile of the airfoil we simply mirror this function horizontally.

# Computing the signed distance function
The signed distance function from the point $p$ to the NACA is the distance from $p$ to the point in the NACA that is closest to it (let's call it $n=[n_x, n_y]^T$).

Since this point is located in the NACA's profile, it follows that 
\\[
\begin{aligned}
n_x &= x \\
n_y &= y_t(x, t)
\end{aligned}
\\]
and since it is the closest to $p$ we know that the distance between the two points should be minimal, with the distance being:
\\[
\operatorname{dist_{n}}=\sqrt{(n_x-p_x)^2+(n_y-p_y)^2}
\\]

Minimizing $\operatorname{dist_{n}}$ is equivalent to minimizing $\operatorname{dist_{n}}^2$, and in order to do that the first thing we need to know is its first derivative: it will be zero when $\operatorname{dist}^2$ is minimal.

The derivative of $\operatorname{dist_{n}}^2$ with respect to $x$ is:
\\[
\frac{d}{dx}\operatorname{dist_{n}}^2=- 2 p_x + 2 x - 2 \left(p_y - yt\right) \frac{d}{d x} y_t
\\]
with:
\\[
\frac{d}{d x}y_t=- 2.074 t x^{3} + 4.2645 t x^{2} - 3.516 t x - 0.63 t + \frac{0.74225 t}{\sqrt{x}}
\\]
Finding the roots of this equation analytically is, unfortunately, either impossible or simply beyond my mathematics skills. We can, however, approximate them numerically by using the *Newton-Rhapson* method.

# Finding roots numerically
## Newton-Rhapson
The [Newton-Rhapson](https://en.wikipedia.org/wiki/Newton%27s_method "Newton-Rhapson") method is an algorithm for approximating roots of real-valued continuously differentiable functions.

The basic idea is that given an initial guess $x_0$, we can approximate a root of a function $f(x)$ by iterating such as:
\\[
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
\\]
## Using N-R to minimize the distance
In our particular case we want to approximate a root of $\frac{d}{dx}\operatorname{dist_{n}}$ so:
\\[
\begin{aligned}
f(x) &= \frac{d}{dx}\operatorname{dist_{n}}^2 \\
f'(x) &= \frac{d^2}{dx^2}\operatorname{dist_{n}}^2
\end{aligned}
\\]
where:
\\[
\frac{d^2}{dx^2}\operatorname{dist_{n}}^2=- 2 \left(p_y - yt\right) \frac{d^{2}}{d x^{2}} y_t + 2 \left(\frac{d}{d x} yt\right)^{2} + 2
\\]
with:
\\[
\frac{d^{2}}{d x^{2}} y_t = \frac{t \left(x^{\frac{3}{2}} \left(- 6.222 x^{2} + 8.529 x - 3.516\right) - 0.371125\right)}{x^{\frac{3}{2}}}
\\]
It should be noted that the roots of the derivative of the distance happen for both extrema of the distance function (both where the distance is minimal and where it is maximal). However, if our initial $x_0$ is close enough to the minimal the Newton-Rhapson method should converge to that (if at all!).

A good initial guess is simply using $x_0=p_x$. This guess is exact for a NACA with thickness $t=0$ (which degenerates the profile into a horizontal segment, and the shortest distance is perpendicular to it) and it gets worse the higher the thickness. It remains a good guess until around $t \sim 0.5$.

After some iterations, we will get a good approximation $x_{\text{approx}}$ such that $\frac{d}{dx}\operatorname{dist_{n}}^2(x_n) \approx 0$, and we can get the corresponding point in the NACA by using the relation:
\\[
n=[x_{\text{approx}}, y_t(x_{\text{approx}}, t)]^T
\\]

## Dealing with x=0
Unfortunately for us, at $=0$ the derivative of the distance function $\frac{d}{dx}\operatorname{dist_{n}}^2$ is not defined due to the $\sqrt{x}$ term in $\frac{d}{dx}{y_t}^2$.

To solve this, we clamp $x$ to be in the $[x_{\epsilon}, 1]$ interval where $x_{\epsilon}$ is a value as close to zero as numerical accuracy allow us. We will call the closest point to $p$ in this "open" NACA $\tilde n$, and the distance to it (that we can compute with Newton-Rhapson as we just saw) $\operatorname{dist_{\tilde n}}$.

Since we just created a gap in the leading edge of the NACA, we should close it. Being as it is a really small gap, we can simply approximate it by placing a segment whose vertices are $(x_\epsilon, y_t(x_{\epsilon}, t))$ and $(x_{\epsilon}, -y_t(x_{\epsilon}, t))$.

The distance from $p$ to this segment can be computed using the signed distance function of a segment (available [here](https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm "here")) and we will call it $\operatorname{dist_{s}}$.

Finally, the distance from the point $p$ to our profile will then be the distance to the (incomplete)NACA or the distance to the gap-filling segment, whichever is smaller, so:
\\[
\operatorname{dist_{n}}=\min({\operatorname{dist_{\tilde n}}, \operatorname{dist_{s}}})
\\]

# Further work / ideas
## High thicknesses
For high thicknesses (higher than around 50%) the Newton-Rhapson method does not converge for the interior of the NACA with our initial guess $x_0=p_x$. This makes sense since the curvature of the leading edge of the NACA deviates the closest point to $p$ to a lower $x$.

In fact, we can actually visualize this by plotting the closest $x$ to every $p=[x,0]^T$ in the [0,1] interval. We can compute this simply by brute force using python. The resulting function looks like this:

![Plot](plot.png "Plot"){: .normal }

with more transparent values representing lower NACA thicknesses.

We can see that for higher thicknesses the distance function is not smooth at all, so it is not surprising that the Newton-Rhapson method has difficulties approximating its root. 

It may still be possible to find a heuristic value for $x_0$ that helps Newton-Rhapson to converge in these situations. However, thicknesses above 40% are already rare in reality, so while the exercise should be interesting, it is probably not that useful.

## Cambered NACA

![NACA8422](NACA8422.png "NACA8422"){: .normal }

There are other, more complex NACA families that are more interesting than the symmetric 00xx series. However, their functions and derivatives are far more involved, and their shapes are so complex that Newton-Rhapson is bound to fail.

By "involved" we mean something such as [this](http://javiermeseguer.com/wp-content/uploads/2021/04/naca-symbolic.html "this").

To overcome this problem we have several alternatives. The easiest one is simply to approximate the NACA by a polygon. This is trivial to do, simply sample the $x$ chord parameter to get $N$ samples and then create segments from $[x_i, \operatorname{naca}(x_i)]^T$ to $[x_i, \operatorname{naca}(x_(i+1)]^T$ for both the upper and lower functions of the desired NACA profile. Then the signed distance function of the NACA can be approximated by that of the polygon.

This is the approximation that we have followed in the following shader:
<iframe width="640" height="360" frameborder="0" src="https://www.shadertoy.com/embed/NdXSz8?gui=true&paused=false&muted=false" allowfullscreen></iframe>

Another choice would be to approximate the NACA by using [B-Splines](https://es.wikipedia.org/wiki/B-spline "B-Splines"). An interesting paper on the topic is ["Universal Airfoil Parametrization Using B-Splines" Dev Rajnarayan, 2018](https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=3150&context=facpub "Universal Airfoil Parametrization Using B-Splines. Dev Rajnarayan, 2018").