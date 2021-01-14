# Complex Plotter

While working through the coursera module [Introduction to Complex Analysis](https://www.coursera.org/learn/complex-analysis/home/welcome) I wanted to develop a robust MATLAB class that allows me to visualize and gain intuition for complex valued functions. 

The comp.plotter class is the result of this ambition.

To instantiate a plotter class, we simply pass a function handle to the constructor. The function handle can be defined or anonymous, there is no difference.

```MATLAB   

my_fun = @(z) z^4 + 10;
p = comp.plotter(my_fun) % Equivalent to p = comp.plotter(@(z) z^2 + 1) 

```
Upon object creation, a series of plots will automatically be generated. Let's examine them.

![](./media/z_4.png)