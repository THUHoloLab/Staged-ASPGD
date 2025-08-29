function MSD = Fx_MSD(varargin)
F3d=varargin{1};
T3d=varargin{2};

MSD=(mean((F3d - T3d).^2,'all'));

end

