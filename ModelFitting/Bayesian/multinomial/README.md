
Multinomial logistic regression in Stan
=======================================

Files included:

-   multinomial\_stan\_comparisons.R: Shows a variety of ways one could write a basic multinomial model with an eye toward speed comparisons.
-   mnl\_1 to 4: stan scripts used in the previous.
-   stan\_basic\_mnl\_conceptual.stan: This is straightforward conceptual stan code for a basic multinomial model. It duplicates the initial example in the manual. As noted there, without a suitable prior it is not identified
-   stan\_basic\_mnl\_altspecific.stan: This extends the standard multinomial model to one that can include both class specific and class constant effects.
