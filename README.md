**Modu** is a mathematical module for doing modular arithmetic in Python. Using operator overloading, it provides an experience similar to a _domain-specific language_ (DSL). Modu is best used as an interactive calculator in Python terminal sessions (REPL) or in Jupyter Notebooks with instantaneous rendering as LaTeX formulas. Modu can be used as well in your Python scripts to ease calculations involving modular arithmetic. The prime target domain is research and education in number theory.

**In a nutshell**, Modu allows you defining and manipulating sets of residue classes for given moduli. These sets are displayed with usual notations of modular arithmetic, using character strings, like `n ≡ {0, ±2, +3} (mod 6)` or using LaTeX formula rendering. New sets can be computed from existing ones using set operations (union, intersection, complement), as well as  arithmetic operations (addition, negation, multiplication, division, exponentiation). Samples of integers can be obtained, displayed as tables (aligning same residues in same columns) and possibly transforming the elements by user-defined functions. Several key results of modular arithmetic can be easily illustrated, like **Fermat's little theorem**, **Euler's totient function** and **Dirichlet's theorem on arithmetic progressions**.

Modu is an open-source module distributed under the **MIT license**.

Here are few examples.
* integers having residue 0 modulo 2 (the even numbers):
```
>>> from modu import O
>>> O % 2
0 (mod 2)
```
* integers multiples of 2 or 3, defined by set union:  
```
>>> O%2 | O%3
{0, ±2, +3} (mod 6)
```
* integers not multiples of 2 or 3, defined by complement of the previous set:
```
>>> ~(O%2 | O%3)
±1 (mod 6)
```
* testing whether 15 and 23 belong to the previous set (note: any prime number greater or equal to 5 belongs to this set!):
```
>>> 15 in ~(O%2 | O%3)
False
>>> 23 in ~(O%2 | O%3)
True
```
* integers multiples of 2 and 3 altogether, defined by set intersection--see Chinese remainder theorem:
```
>>> O%2 & O%3
0 (mod 6)
```
* checking whether two sets are equivalent or not:
```
>>> O%2 | O%3 == (0, 2, 3, 4) + O%6
True
>>> O%2 & O%3 == O%6
True
>>> O%2 | O%3 != O%6
True
```
Modu provides more functions for normalizing the representation, changing the display format or providing samples using table layout.

All those features are demonstrated in an interactive tutorial, the **"_Modutorial_" Jupyter Notebook** (modutorial.ipynb). This can be opened directly online by clicking here (binder):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/piedenis/Modu/HEAD?urlpath=%2Fdoc%2Ftree%2Fsrc%2Fmodutorial.ipynb)

(be patient: you have to wait a bit...)
