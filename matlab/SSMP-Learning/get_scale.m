function [ dataDecayNorm , lowDimCtrlNorm, obj ] = get_scale( data, dataDecay, lowDimCtrl , params )
    %scale: Scale sim/exp data to be in range [-1 , 1]
    %    Also creates scaleup/scaledown matrices and saves as params
    %    data - struct containing fields t , y , u (at least)
    %    data_scaled - struct containing t , y , u , x (optional)   
    
    % get min/max values in each dimension
    y_min = min( data.y );
    y_max = max( data.y );
    u_min = min( data.u );
    u_max = max( data.u );
    
    % calculate centers of range
    y_dc = ( y_max + y_min ) ./ 2;
    u_dc = ( u_max + u_min ) ./ 2;
    
    % calculate scaling factors
    scale_y = ( y_max - y_min ) ./ 2;
    scale_u = ( u_max - u_min ) ./ 2;
    
    % shift and scale the data
    data_scaled = struct;    % initialize
    data_scaled.t = data.t;  % time is not scaled
    data_scaled.y = ( data.y - y_dc ) ./ scale_y;
    data_scaled.u = ( data.u - u_dc ) ./ scale_u;
    
    % save scaling functions
    y = sym( 'y' , [ 1 , params.sizeX ] );
    y_scaledown = ( y - y_dc ) ./ scale_y;
    obj.scaledown.y = matlabFunction( y_scaledown , 'Vars' , {y} );
    
    y_scaleup = ( y .* scale_y ) + y_dc;
    obj.scaleup.y = matlabFunction( y_scaleup , 'Vars' , {y} );

    u = sym( 'u' , [ 1 , params.sizeU ] );
    u_scaledown = ( u - u_dc ) ./ scale_u;
    obj.scaledown.u = matlabFunction( u_scaledown , 'Vars' , {u} );

    u_scaleup = ( u .* scale_u ) + u_dc;
    obj.scaleup.u = matlabFunction( u_scaleup , 'Vars' , {u} );
    
    % save scaling factors
    obj.params.scale.y_factor = scale_y;
    obj.params.scale.y_offset = y_dc;
    obj.params.scale.u_factor = scale_u;
    obj.params.scale.u_offset = u_dc;

    dataDecayNorm = cell(size(dataDecay,2),2);
    for i = 1:size(dataDecay,2)
        dataDecayNorm{i, 1} = (dataDecay{1, i}.t)';
        dataDecayNorm{i, 2} = (obj.scaledown.y(dataDecay{1, i}.y))';
    end
    
    lowDimCtrlNorm = cell(size(lowDimCtrl,2),2);
    for i = 1:size(lowDimCtrl,2)
        lowDimCtrlNorm{i, 1} = (lowDimCtrl{1, i}.t)';
        lowDimCtrlNorm{i, 2} = (obj.scaledown.y(lowDimCtrl{1, i}.y))';
        lowDimCtrlNorm{i, 3} = (obj.scaledown.u(lowDimCtrl{1, i}.u))';
    end

end