While revising my concepts of pricing various financial models, I came across the fact that internet doesn't distinguish between clean/dirty prices while computing Internal rate of return for a bond.

After some experimentation, I finally figured out the formula to compute the same. The basic relation is:

- Market quotes the price as `clean price`.
- Clean price to dirty price can be computed as: Dirty Price = Clean Price + Coupon * fraction of time passed between coupons. (Final variable doesn't depend on market quoted yield).
- Fraction of time is computed via `Day Count Convention`.
- IRR is computed via clean price.
