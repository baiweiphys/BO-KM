``BO-KM`` (available at the URL: https://github.com/baiweiphys/BO-KM/) is a new tool developed by Wei Bai and Huasheng Xie. Its primary objective is to efficiently calculate all solutions for the dispersion relation of obliquely propagating waves in magnetized hot plasmas with multiple species. This code provides extensive support for a wide range of distribution functions, including anisotropic drift kappa-Maxwellian distributions, anisotropic drift bi-Maxwellian distributions, and their combinations. Additionally, it can also handle bi-kappa distributions for perpendicular propagation.

Wei Bai (Email: baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn)
College of Electrical and Power Engineering, Taiyuan University of Technology, Taiyuan 030024, China

Huasheng Xie (Email: huashengxie@gmail.com, xiehuasheng@enn.cn)
ENN Science and Technology Development Co., Ltd., Langfang 065001, China,
Hebei Key Laboratory of Compact Fusion, Langfang 065001, China

2023-11-16

## BO-KM

``BO-KM`` is a comprehensive code that has the capability to simultaneously solve for all significant roots of the dispersion relation for obliquely propagating waves in magnetized multi-species plasma. The code supports various distributions, including anisotropic drift kappa-Maxwellian distributions, anisotropic drift bi-Maxwellian distributions, and a combination of both (as well as bi-kappa distributions for perpendicular propagation). The ``BO-KM`` solver provides a powerful and efficient tool for analyzing the waves and instabilities in magnetized multi-species plasma with anisotropic drift kappa-Maxwellian, bi-Maxwellian distributions, or a mixture of both.

Input: The background magnetic field is denoted as $B_0$ (in Tesla), the angle of oblique propagation is represented as deg (in degrees), $J$ is a non-negative integer for $J$-pole expansion, $N$ is a non-negative integer for Bessel functions, and nk indicates the number of grid points.
The ``bokm.in`` file requires the following input parameters: $q_s$ (charge of species in elementary charge units), $m_s$ (mass of species in terms of electron mass), $n_s$ (number density of species in cubic meters), $T_{\parallel s}$ (parallel temperature of species in electron volts), $T_{\perp s}$ (perpendicular temperature of species in electron volts), $u_{s_0}/c$ (species drift velocity normalized by the speed of light), $\kappa$ (value of kappa for the distribution function of the $s$-th particle species), and $\kappa_{s,th}$ (threshold value for kappa).

We present seven benchmark examples that serve as representative cases, including the Multi-species bump-on-tail mode, EMEC waves, Bernstein waves, Whistler waves, Mirror mode, and EMIC waves.

At present, only Matlab version of ``BO-KM`` is available, which has been tested in MacOS and Windows 11 under Matlab 2023b.

## Files structures:

- ``BO-KM/RUN/Case_XXX/``    %  Benchmark  cases 
	- ``BO-KM/RUN/Case_XXX/main_bokm.m``    % Driver routine of BO-KM. Execute this code to initiate the simulation.
	- ``BO-KM/RUN/Case_XXX/output`` % For storing the output data.
		- ``selected_plot``  % A subroutine is developed by Xie (Xie 2019, CPC) to enhance the visualization of 2D or 3D plots, enabling the selection of a specific dispersion surface or curve.
- BO-KM/kernels                        %  Kernel file for ``BO-KM``
	- ``solver_km.m``                   % Solver for kappa-Maxwellian plasmas
	 - ``solver_maxwell.m``        % Solver for bi-Maxwellian plasmas
	 - ``solver_mixed.m``            % Solver for mixed distributions in both kappa-Maxwellian and bi-Maxwellian plasmas.
	 - ``func_Jpole.m``                % Function for J-pole expansion (from Xie2016, PST)
	 - ``getIndexOfBlkMatrix_km.m``  % To obtain the index of sub-block Matrix in KM total Matrix.
	 -  ``getIndexOfBlkMatrix_maxwell.m``    % To obtain the index of sub-block Matrix in the BM total Matrix
	 - ``getIndexOfBlkMatrix_mixed.m``          % To get the index of sub-block Matrix in the total Matrix mixed with BK and KM.
	 - ``maxwell_bsnj.m``, ``maxwell_bx10.m``, ....  % The coefficients for the oblique bi-Maxwellian plasma model.
	 - ``km_bsl.m``, ``km_bx10.m``,...  % The coefficients for the oblique kappa-Maxwellian plasma model.
	 - ``M_maxwell.m``   % Equivalent matrix for bi-Maxwellian plasma model.
	 - ``Ml_km.m``,``Mlm1_km.m``, ``Mlm2_km.m``, ``Mlp1_km.m``  % The four equivalent matrices for kappa-Maxwellian plasma model.
	 - ``Ml_mixed_km.m``,``Mlm1_mixed_km.m``, ``Mlm2_mixed_km.m``, ``Mlp1_mixed_km.m``, ``M_mixed_maxwell.m``  % The five equivalent matrices for the plasma model mixed with kappa-Maxwellian and bi-Maxwellian distributions.
	 - ``getLen_Ml.m``, ``getLen_Mlm1.m``, ``getLen_Mlm2.m``, ``getLen_Mlp1.m``  % Functions to obtain the size of the matrices Ml_km, Mlm1_km, Mlm2_km, and Mlp1_km.
	 - ``getPlasmaPrameters.m``, % To obtain the plasma parameters.
- BO-KM/PerpSolver4BK  % Sub-module for perpendicular propagation in bi-kappa plasmas.
	- ``biKappa_bx1sn.m``,``biKappa_bx2sn.m``,``biKappa_by1sn.m``,``biKappa_by2sn.m``,``biKappa_bz3.m``,``biKappa_bz3sn.m``  % Calculate the coefficients for the perpendicular plasma waves with anisotropic drift bi-kappa distribution.
	- ``getIndexOfBlkMatrix_bikappa.m`` % To obtain the index of a subblock matrix within a matrix for perpendicular propagation in bi-kappa plasmas.
	- ``M_bikappa.m``  % Compute the matrix M_bikappa for the perpendicular propagation plasma wave model with a bi-kappa distribution.
	- ``S0snFunc.m``,``S1snFunc.m``,``S3snFunc.m``,``S5snFunc.m``  % Compute the integral of S0sn, S1sn, S3sn, and S5sn, which are the expression obtained from Ref. (Summers1994, PoP) .
	- ``solver_bikappa.m``   %To calculate the roots for the perpendicular propagation plasma wave model with a bi-kappa distribution, given the value of kx and kz.
- BO-KM/toolBox/``createDateFile.m`` % Create a new data folder if it doesn't exist.

## How to run the code:

1. 1. Set input parameters in the file ``BO-KM/RUN/case_XXXX/bokm.i``.
2. Set the parameters of ``N``, ``J``, ``deg`` (propagation angles, in degrees), and ``B0`` (background magnetic field in the z direction) in the file ``BO-KM/RUN/case_XXXX/main_bokm.in``.
3. After seting the above parameters in ``main_bokm.m``, run it. You can see some plots of the results. And all results will be stored in output directory. You can modify the plot section in the file  for different plots, or writing your own plot files.
4. After setting the above parameters in the file ``main_bokm.m``, run it. You can see some plots of the results. Additionally, all the output will be saved in the ``output`` directory. If you wish to customize the plots, you can modify the corresponding section in the file or create your own plot files.

If you encounter any issues, kindly get in touch with us.

Wei Bai, baiweiphys@gmail.com, baiwei12@mail.ustc.edu.cn, TYUT

Huasheng Xie,  huashengxie@gmail.com, xiehuasheng@enn.cn, ENN

2023-11-16 


