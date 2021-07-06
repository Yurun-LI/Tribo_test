# the simulation of tribometer test
------

## 1.	Introduction

This study simulates the frictional behavior of a tribological conjunction lubricated with nanometer-thick perfluoropolyether (PFPE) lubricant. In the simulation, a magnetic disk, composing of glass substrate, magnetic layer and 3 nm thick diamond-like carbon overcoat, is rotated against a BK-7 glass pin with spherical end cap of 7.5 mm radius. The initial supplied lubricant, hini along this sliding contact consists of 1.7 nm thick non-polar Z03 type PFPE lubricant. The simulated tribological conjunction is illustrated in following figure.

![model.jpg](https://github.com/Yurun-LI/Tribo_test/blob/master/Images/model.jpg?raw=true)

## 2.	Mathematic Model

> the model is built based on Elrod's Cavitation Algorithm, and the assumptions and conditions is as following $\mathrm{Eqn. 1}$:

### 2.1	2D Reynold's Equation for point contact

$$
{\partial\over{x}}[{{\rho}h^{3}\over{\eta}}\cdot{\partial{p}\over{\partial{x}}}]+{\partial\over{y}}[{{\rho}h^{3}\over{\eta}}\cdot{\partial{p}\over{\partial{y}}}] = 12{\{{\partial\over\partial{x}}[{{\rho}{h}{(u_{av})}}]\}} 
$$

### 2.2	Pressure Equation

the pressure equation as derived by Elrod is as shown in $\mathrm{Eqn. 2}$
$$
p = g{\beta}{\ln}{\theta}+p_{c}, where {\space}{\space}{\space}{\space}\theta\ne0
$$

### 2.3	Density Equation

Erode relates the frictional film content, $\theta$ to the ratio of the density at full film region and cativation region. the relationship can be expressed as $\mathrm{Eqn. 3}$
$$
{\theta} = {{\rho}\over{\rho_c}}\\{{\rho_c}={{\rho}\over{\theta}}}
$$
