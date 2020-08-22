# Michaelis-Menten Flask Calculator

Determine <img src="https://render.githubusercontent.com/render/math?math=k_\mathrm{cat}"> (turnover number measuring the rate of the enzyme) and <img src="https://render.githubusercontent.com/render/math?math=K_\mathrm{M}"> (Michaelis constant measuring subtrate binding affinity) from enzyme activity assay. 
Performs hyperbolic curve fitting with <img src="https://render.githubusercontent.com/render/math?math=V_0"> (initial velocity data) at varying substrate concentrations to the Michaelis-Menten equation.  

<img src=
"https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AV_0+%26%3D+%5Cfrac%7BV_%5Cmathrm%7Bmax%7D+%5Ccdot+%5B%5Cmathrm%7BS%7D%5D%7D+%7BK_%5Cmathrm%7BM%7D+%2B+%5B%5Cmathrm%7BS%7D%5D%7D+%5C%5C%0A+%26%3D+%5Cfrac%7Bk_%7Bcat%7D+%5Ccdot+%5B%5Cmathrm%7BE%7D%5D+%5Ccdot+%5B%5Cmathrm%7BS%7D%5D%7D%7BK_%5Cmathrm%7BM%7D+%2B+%5B%5Cmathrm%7BS%7D%5D%7D%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
V_0 &= \frac{V_\mathrm{max} \cdot [\mathrm{S}]} {K_\mathrm{M} + [\mathrm{S}]} \\
 &= \frac{k_{cat} \cdot [\mathrm{E}] \cdot [\mathrm{S}]}{K_\mathrm{M} + [\mathrm{S}]}
\end{align*}
">

When <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cbegin%7Balign%2A%7D%0AK_%7B%5Cmathrm%7BM%7D%7D+%5Cgg+%5B%5Cmathrm%7BS%7D%5D%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
K_{\mathrm{M}} \gg [\mathrm{S}]
\end{align*}
">, we utilize the simplified form of the the Michaelis-Menten equation and solve for <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cbegin%7Balign%2A%7D%0Ak_%7B%5Cmathrm%7Bcat%7D%7D%2FK_%7B%5Cmathrm%7BM%7D%7D+%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
k_{\mathrm{cat}}/K_{\mathrm{M}} 
\end{align*}
"> (catalytic efficiency) instead of <img src="https://render.githubusercontent.com/render/math?math=k_\mathrm{cat}">  and <img src="https://render.githubusercontent.com/render/math?math=K_\mathrm{M}"> separately.

<img src=
"https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AV_0+%3D+%5Cfrac%7Bk_%7B%5Cmathrm%7Bcat%7D%7D+%5Ccdot+%5B%5Cmathrm%7BE%7D%5D+%5Ccdot+%5B%5Cmathrm%7BS%7D%5D+%7D+%7BK_%7B%5Cmathrm%7BM%7D%7D%7D%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
V_0 = \frac{k_{\mathrm{cat}} \cdot [\mathrm{E}] \cdot [\mathrm{S}] } {K_{\mathrm{M}}}
\end{align*}
">

## Installation
```
git clone https://github.com/eking2/enz_kinetics_flask.git
cd enz_kinetics_flask
pipenv shell
pipenv install
python app.py
```

## Example

### Input

Enter the listed reaction parameters. 
The slopes are entered as space delimited format with the first column as substrate concentrations in mM. 
The next columns are the absorbance slopes in mAbs/min, this can be directly copied from Excel. 
A set of toy data is provided in the `test_input` folder. 
Single trial or multiple trials can be fit at once as long as the trials are separated into different columns.

![input](assets/inputs.png)

### Output

The output is a Michaelis-Menten plot with the calculated kinetic parameters. The complete raw input and output data can be downloaded from the link. 

![input](assets/output.png)
