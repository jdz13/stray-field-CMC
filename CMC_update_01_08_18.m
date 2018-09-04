% Coherent code for all of the outputs. Should allow for quick changes in variable that will output all results, with no additional work.
%
% made by JDZ 31/10/2017


% One method utilises the Akoun method in 3D - otherwise found by Janssen
% in their PhD thesis. Nothing special to note.
%
% One method utilises the Dipole approximation to find the stray field at
% a distance B. The method is taken from Fundamentals of Magnetism by Mario
% Reis - chapter 1. Nothing special to note. 
%
% One method uses the scalar potential method to obtain the stray field at
% a distance. The method is taken from Engel-Herbert & Hesjedal (2005)
% Care should be taken, as this uses equations where the magnetisation is
% in y (in plane) rather than z (out of plane). Switch dimensions. 
%
% One method uses the far field equation
%
% Lastly, values are pulled from MuMax calculations and extrapolated


% Integration space defined by (xa,ya,za)
% Magnetic region defined by (xb,yb,zb)
% Muvalue is from experimental VSM results from a 1nm CoFeB sample
% 
%
% Yeilds B at a distance. 
%
% JDZ Sep 2017
%

clear

%Define constants

mu0 = 4* pi()* 10^-7; % [H/m]

%Define the number of divisions and linspace variables

div = 50; % total number of divisions so, div/2 on each side (and one in centre)

%Define the magnetic area

dist = 10^-3; % distance in [m] in both positive and negative directions. 
xb = dist/div/2; % [m] % div by 2 for factor 8 later? 
yb = dist/div/2; % [m]
zb = dist/div/2; % [m]

%Define the minimum space to be integrated over, what we would normally
%refer to the integration space. 

xa = dist; % [m]
ya = dist; % [m]
za = dist; % [m]

    %Dipole expression breaks dipoles up, so seperate division needed
    div3d = div; 
    divmag = 1;
    linsp_dipole = linspace(-za,za,div3d+1); % this has to be above the loop 

%Work out the division spacing, and in turn then make the grid. 

divspacing = [2*xa/div,2*ya/div,2*za/div]; % just good to know

Xm = linspace(-xa,xa,div+1);
Ym = linspace(-ya,ya,div+1);
Zm = linspace(-za,za,div+1);
linspAkoun = Zm;

[X,Y,Z] = meshgrid(Xm,Ym,Zm);

% Calculating M 
M = 10^6; %  
OurA = 8*xb*yb*zb;
Muvalue = M*OurA; 

    %For Magnetostatic_H this is enough

    %For Akoun, we need to define Br
    Br = M;
    const = Br/(4*pi());
    
    %For Dipole Expression, we need to ensure Muvalue is in Z direction,
    %descretised for each dipole
    NormalisedMu = (Muvalue/(divmag+1)^3);
    mu = [0,0,NormalisedMu];
        
%Define the H variables (empty), for code efficiency

HxAkoun = zeros(size(X,1), size(X,2), size(X,3));
HyAkoun = HxAkoun;
HzAkoun = HxAkoun;

HxMagH = zeros(size(X,1), size(X,2), size(X,3));
HyMagH = HxMagH;
HzMagH = HxMagH;

        %For the dipole expression
        Xmdipole = linspace(-xa,xa,div3d+1);
        Ymdipole = linspace(-ya,ya,div3d+1);
        Zmdipole = linspace(-za,za,div3d+1);

        [Xdipole,Ydipole,Zdidpole] = meshgrid(Xmdipole,Ymdipole,Zmdipole);

        %Define the H variables (empty), for code efficiency

        Hxdipole = zeros(size(Xdipole,1), size(Xdipole,2), size(Xdipole,3));
        Hydipole = Hxdipole;
        Hzdipole = Hxdipole;

%%Start to look at getting the equations solved for in all space

%AKOUN METHOD
for a = linspace(-xa,xa,div+1)
    for b = linspace(-ya,ya,div+1)
        for c = linspace(-za,za,div+1)
          
            
            for k = 0:1
                for l = 0:1
                    for m = 0:1
                        
                        S = a - ((-1)^k)*xb;
                        T = b - ((-1)^l)*yb;
                        U = c - ((-1)^m)*zb;
                        R = sqrt(S^2+T^2+U^2);
                 
                        dHzAkoun = ((-1)^(k+l+m))* atan((S*T)/(U*R));
                        
                        unit1 = round((a+xa)*(div/2/xa))+1;
                        unit2 = round((b+ya)*(div/2/ya))+1;
                        unit3 = round((c+za)*(div/2/za))+1;
                        
                        HzAkoun (unit1,unit2,unit3) = HzAkoun(unit1,unit2,unit3) + (const*dHzAkoun);
                     
                    end
                end 
            end
          
            %MagnetostaticH method
            for k = 1:2
                for l = 1:2
                    for m = 1:2
                        
                        Xfac = a + ((-1)^k)*xb;
                        Yfac = c + ((-1)^l)*zb;
                        Zfac = b + ((-1)^m)*yb;
                        SqrtXYZ = sqrt(Xfac^2+Yfac^2+Zfac^2);
                        
                        dHxMagH = (-1)^(k+l+m)*log(Zfac+SqrtXYZ);
                      
                        dHzAMagH = (Yfac*Xfac)/((abs(Yfac))*(abs(Xfac)));
                        dHzBMagH = atan((dot(abs(Xfac),Zfac))/(dot(abs(Yfac),SqrtXYZ)));
                        dHzMagH = ((-1)^(k+l+m))*dHzAMagH*dHzBMagH;
                        
                        dHyMagH = (-1)^(k+l+m)*log(Xfac+SqrtXYZ);
                        
                        unit1 = round((a+xa)*(div/2/xa))+1;
                        unit2 = round((c+za)*(div/2/za))+1;
                        unit3 = round((b+ya)*(div/2/ya))+1;
                        
                        HxMagH (unit1,unit2,unit3) = HxMagH(unit1,unit2,unit3) + 1/4/pi()*dHxMagH*M;
                        HzMagH (unit1,unit2,unit3) = HzMagH(unit1,unit2,unit3) - 1/4/pi()*dHzMagH*M;
                        HyMagH (unit1,unit2,unit3) = HyMagH(unit1,unit2,unit3) + 1/4/pi()*dHyMagH*M;
                     
                    end
                end 
            end
            
        end 
    end 
end 

% Manipulate the data to compare with other methods

% AKOUN METHOD

BzAkoun = HzAkoun.*(4*pi()*10^-7);
middleakoun = zeros(1,size(BzAkoun,3));

BxMagH = HxMagH.*mu0;
ByMagH = HyMagH.*mu0;
BzMagH = HzMagH.*mu0;
middleMagH = zeros(1,size(BzMagH,3));

for o = 1:size(BzAkoun,3)  
    go = BzAkoun(:,:,o);
    go = reshape(go,size(BzAkoun,1),size(BzAkoun,2));
  
    middleakoun(o) = go(round((div+1)/2),round((div+1)/2));
    
    go = BzMagH(:,o,:);
    go = reshape(go,size(BzMagH,1),size(BzMagH,2));
  
   middleMagH(o) = go(round((div+1)/2),round((div+1)/2));

end 

for a = linspace(-xa,xa,div3d+1)
    for b = linspace(-ya,ya,div3d+1)
        for c = linspace(-za,za,div3d+1) 
            
            for d = linspace(-xb, xb, divmag+1)
                for e = linspace(-yb, yb, divmag+1)
                    for f = linspace(-zb, zb, divmag+1)
                    
                        r = [(a-d), (b-e),(c-f)];
                        rscalar = sqrt((a-d)^2 + (b-e)^2 +(c-f)^2);
                        
                        Hdipole = ((3*dot(mu,r)*r/(rscalar^5))-(mu/(rscalar^3)))/4/pi();
                        
                        unit1 = round((a+xa)*(div3d/2/xa))+1;
                        unit2 = round((b+ya)*(div3d/2/ya))+1;
                        unit3 = round((c+za)*(div3d/2/za))+1;
                        
                        Hxdipole (unit1,unit2,unit3) = Hxdipole(unit1,unit2,unit3) + Hdipole(1);
                        Hydipole (unit1,unit2,unit3) = Hydipole(unit1,unit2,unit3) + Hdipole(2);
                        Hzdipole (unit1,unit2,unit3) = Hzdipole(unit1,unit2,unit3) + Hdipole(3);
                     
                    end
                end 
            end
          
            
        end 
    end 
end 

% need some form of graphical output, to be examined 

Bxdipole = Hxdipole.*(4*pi()*10^-7);
Bydipole = Hydipole.*(4*pi()*10^-7);
Bzdipole = Hzdipole.*(4*pi()*10^-7);

% Now work out all of the points in the centre for the plot. 

middledipole = zeros(1,size(Bzdipole,3));

for o = 1:size(Bzdipole,3)  
    go = Bzdipole(:,:,o);
    go = reshape(go,size(Bzdipole,1),size(Bzdipole,2));
    
    middledipole(o) = go(round((div3d+1)/2),round((div3d+1)/2));
    
end 

Far_Field = 2*mu0/4/pi()*(Zm.^-3)*Muvalue;

Full = sqrt(Bzdipole.^2+Bydipole.^2+Bxdipole.^2);
%%
figure(1)
clf
semilogy(linsp_dipole, middledipole, 'ko', Zm, Far_Field,'r-', Zm,middleakoun,'bx', Zm, middleMagH,'g--')
legend('Dipole Expression', 'Far Field Expression', 'Akoun 3d', 'Magnetostatic H')
title ('Stray field at a distance')
xlabel ('Distance in Z [m]')
ylabel('Stray field - B_Z [T]')
axis([0 0.001 10^(-9) 10^(-2) ])

%%
figure(2)
clf
slice(X,Y,Z,Full, [],0,[])
caxis([0,10^-4])