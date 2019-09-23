function mu = SW_Viscosity(T,uT,S,uS)
    % SW_Viscosity    Dynamic viscosity of seawater
    %=========================================================================
    % USAGE:  mu = SW_Viscosity(T,uT,S,uS)
    %
    % DESCRIPTION:
    %   Dynamic viscosity of seawater at atmospheric pressure (0.1 MPa) using
    %   Eq. (22) given in [1] which best fit the data of [2], [3] and [4].
    %   The pure water viscosity equation is a best fit to the data of [5].
    %   Values at temperature higher than the normal boiling temperature
    %   are calculated at the saturation pressure.
    %
    % INPUT:
    %   T  = temperature
    %   uT = temperature unit
    %        'C'  : [degree Celsius] (ITS-90)
    %        'K'  : [Kelvin]
    %        'F'  : [degree Fahrenheit]
    %        'R'  : [Rankine]
    %   S  = salinity
    %   uS = salinity unit
    %        'ppt': [g/kg]  (reference-composition salinity)
    %        'ppm': [mg/kg] (in parts per million)
    %        'w'  : [kg/kg] (mass fraction)
    %        '%'  : [kg/kg] (in parts per hundred)
    %
    %   Note: T and S must have the same dimensions
    %
    % OUTPUT:
    %   mu = dynamic viscosity [kg/m-s]
    %
    %   Note: mu will have the same dimensions as T and S
    %
    % VALIDITY: 0 < T < 180 C and 0 < S < 150 g/kg;
    %
    % ACCURACY: 1.5%
    %
    % REVISION HISTORY:
    %   2009-12-18: Mostafa H. Sharqawy (mhamed@mit.edu), MIT
    %               - Initial version
    %   2012-06-06: Karan H. Mistry (mistry@mit.edu), MIT
    %               - Allow T,S input in various units
    %               - Allow T,S to be matrices of any size
    %
    % DISCLAIMER:
    %   This software is provided "as is" without warranty of any kind.
    %   See the file sw_copy.m for conditions of use and licence.
    %
    % REFERENCES:
    %   [1] M. H. Sharqawy, J. H. Lienhard V, and S. M. Zubair, Desalination
    %       and Water Treatment, 16, 354-380, 2010. (http://web.mit.edu/seawater/)
    %   [2] B. M. Fabuss, A. Korosi, and D. F. Othmer, J., Chem. Eng. Data 14(2), 192, 1969.
    %   [3] J. D. Isdale, C. M. Spence, and J. S. Tudhope, Desalination, 10(4), 319 - 328, 1972
    %   [4] F. J. Millero, The Sea, Vol. 5, 3 – 80, John Wiley, New York, 1974
    %   [5] IAPWS release on the viscosity of ordinary water substance 2008
    %=========================================================================

    %% CHECK INPUT ARGUMENTS

    % CHECK THAT S&T HAVE SAME SHAPE
    if ~isequal(size(S),size(T))
        error('check_stp: S & T must have same dimensions');
    end

    % CONVERT TEMPERATURE INPUT TO °C
    switch lower(uT)
        case 'c'
        case 'k'
            T = T - 273.15;
        case 'f'
            T = 5/9*(T-32);
        case 'r'
            T = 5/9*(T-491.67);
        otherwise
            error('Not a recognized temperature unit. Please use ''C'', ''K'', ''F'', or ''R''');
    end

    % CONVERT SALINITY TO PPT
    switch lower(uS)
        case 'ppt'
        case 'ppm'
            S = S/1000;
        case 'w'
            S = S*1000;
        case '%'
            S = S*10;
        otherwise
            error('Not a recognized salinity unit. Please use ''ppt'', ''ppm'', ''w'', or ''%''');
    end

    % CHECK THAT S & T ARE WITHIN THE FUNCTION RANGE
    if ~isequal((T<0)+(T>180),zeros(size(T)))
        warning('Temperature is out of range for Viscosity function 0<T<180 C');
    end

    if ~isequal((S<0)+(S>150),zeros(size(S)))
        warning('Salinity is out of range for Viscosity function 0<S<150 g/kg');
    end

    %% BEGIN

    S = S/1000;

    a = [
        1.5700386464E-01
        6.4992620050E+01
       -9.1296496657E+01
        4.2844324477E-05
        1.5409136040E+00
        1.9981117208E-02
       -9.5203865864E-05
        7.9739318223E+00
       -7.5614568881E-02
        4.7237011074E-04
    ];

    mu_w = a(4) + 1./(a(1)*(T+a(2)).^2+a(3));


    A  = a(5) + a(6) * T + a(7) * T.^2;
    B  = a(8) + a(9) * T + a(10)* T.^2;
    mu = mu_w.*(1 + A.*S + B.*S.^2);

end
