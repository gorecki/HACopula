# The HACopula Toolbox
## Hierarchical Archimedean copulas for MATLAB and Octave

## Description:
The HACopula toolbox extends the copula modeling provided by MATLAB to modeling with hierarchical Archimedean copulas, which allows for non-elliptical distributions in arbitrary dimensions enabling for asymmetries in the tails. This includes their representation as MATLAB objects, evaluation, sampling, estimation and goodness-of-fit testing, as well as tools for their visual representation or computation of corresponding Kendall's tau and tail dependence coefficients. The toolbox is also compatible with Octave, where no support for copulas in more than two dimensions is currently provided.

## Manual:
Available as PDF, see MANUAL.pdf.

## Installation:
It is enough to unpack the files to a selected folder and to add this folder with its subfolders to the MATLAB (Octave) path. For the full functionality in MATLAB, the toolbox requires Statistics and Machine Learning Toolbox and Symbolic Math Toolbox. In Octave, it suffices to install and to load the 'statistics' package from https://octave.sourceforge.io/statistics/index.html and the symbolic package OctSymPy (https://github.com/cbm755/octsympy).

## Examples:
Three examples covering most of the functionality provided by the toolbox are placed in the folder Demos, namely quickex.m, elaboratedex.m and highdimex.m.

## Licence:    
HACopula is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HACopula is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HACopula.  If not, see <http://www.gnu.org/licenses/>


Copyright 2018 Jan Górecki (gorecki@opf.slu.cz)

(Please, report all suggestions, comments, questions or bugs to the above-mentioned e-mail)


