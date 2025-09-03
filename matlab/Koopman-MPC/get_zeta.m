function [ zeta ] = get_zeta( params , data_in )
%get_zeta: Adds a zeta field to a test data struct
%   data_in - struct with t , x , y , u fields
%   zeta - [ y , yd1 , yd2 , ... , ud1 , ud2 , ... ]
sizey = 2*params.n + params.m;
zeta = zeros(1,sizey);
% add the zeta field
for i = params.nd + 1 : size( data_in.y , 1 )
    ind = i - params.nd;    % current timestep index
    y = data_in.y( i , : );
    u = data_in.u( i , : );
    ydel = zeros( 1 , params.nd * params.n );
    udel = zeros( 1 , params.nd * params.m );
    for j = 1 : params.nd
        fillrange_y = params.n * (j - 1) + 1 : params.n * j;
        fillrange_u = params.m * (j - 1) + 1 : params.m * j;
        ydel(1 , fillrange_y) = data_in.y( i - j , : );
        udel(1 , fillrange_u) = data_in.u( i - j , : );
    end
    zetak = [ y , ydel , udel ];
%                 if obj.liftinput == 1     % include input in zeta
%                     zetak = [ zetak , u ];
%                 end
    zeta( ind , : ) = zetak;
    % uzeta( ind , : ) = data_in.u( i , : );    % current timestep with zeta (input starting at current timestep)
end
end