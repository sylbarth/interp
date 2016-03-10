# Interpolation Library

Most of the times, in engineering, science, finance and economic analysis, data is given at discrete points `(x1, ... ,xn)` and `(y1, ... yn)`, and obtained by sampling or experimentation. These datasets represent the values of a function for a limited number of values, but how then does one find the value of y at any other value of x ?

A continuous function `f(x)` may be used to represent the missing data values. An interpolation is therefore a method of constructing these new data points. Then one can find the value of y at any other value of x.

This C++ Library provides three interpolation methods: Lagrande, Divided Differences and Cubic Splines.

