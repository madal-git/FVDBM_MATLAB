function BODY_FORCE = body_force (Gravity, Density, Temperature, Temperature_ref, FTH, FBF)

% Gravity is the gravity force, a column vector
% Density is the fluid density, a row vector
% Temperature is the fluid temperature, a row vector
% Temperature_ref is the reference temperature, a single value
% FTH is the flag for thermal model. FTH=0---isothermal; FTH=1---thermal
% FBF is the flag for different model of body force




if FTH==0% Body force, isothermal
    BODY_FORCE=Gravity*Density;
elseif FTH==1 % Gravitational force, thermal
    beta=1/Temperature_ref;
    %         F_source=beta*(gravity*(Density.*(Temperature-T_ref)));
    BODY_FORCE=-beta*(Gravity*(Density.*(Temperature-Temperature_ref)));
else
    error('Wrong flag for thermal model!');
end