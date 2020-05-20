%CRU2GDL
%this script loads the gridded, monthly weather data of the Climate
%Research Unit of the University of East Anglia and computes the tridecadal
%average temperature for the subnational units of the Global Data Lab
%
%Data source 1
%http://www.cru.uea.ac.uk/cru/data/hrg/
%Make sure you unzip the .dat file
%
%The original data is on a 0.5x0.5 degree grid. The data is organized by
%grid, month, year. That is, the first block of 360x720 data points are for
%January 1900, the second block for February 1900, and so on all the way
%to December 2019.
%
%This procedure works for the temperature data but it can easily be
%adjusted for other climate parameters: Just replace the input file.
%
%The output B is the average annual temperature, gridded, for 1990-2010.
%This can be easily changed too, just watch your indices.
%
%Data source 2
%https://globaldatalab.org/shdi/shapefiles/
%That shapefile comes with the Subnational Human Development Index
%More data for the same regions at https://globaldatalab.org/shdi/gnic/?levels=1%2B4&interpolation=0&extrapolation=0&nearest_real=0
%
%The output T is the average annual temperature for each administrative area.
%
%The script relies on the Matlab Mapping Toolbox
%
%Richard S.J. Tol, 20 May 2020

clear all

NLong = 720;
NLat = 360;
NMonth = 12;
NYear = 30;
Res = 2; %0.5x0.5 degrees makes 2 grid cells per degree

load cru_ts4.04.2001.2010.tmp.dat
load cru_ts4.04.1991.2000.tmp.dat
load cru_ts4.04.1981.1990.tmp.dat

A80=cru_ts4_04_1981_1990_tmp;
A90=cru_ts4_04_1991_2000_tmp;
A00=cru_ts4_04_2001_2010_tmp;
clear cru*
A80(A80==-999)=NaN;
A90(A90==-999)=NaN;
A00(A00==-999)=NaN;

A = A80 + A90 + A00; %we compute the average, so adding is fine

%%
B = zeros(NLat, NLong);
for y = 1:(NYear/3) %we added the three decades
    for m = 1:NMonth
        for l = 1:NLat
            Aind = l+(m-1)*NLat+(y-1)*NLat*NMonth;
            B(l,:) = B(l,:) + A(Aind,:)/NMonth/NYear/10; %divide by 10 to get celsius
        end
    end
end

clear A*
clear y m l NMonth NYear

%%
Long = -180;
for i = 2:NLong,
    Long(i) = Long(i-1) + 0.5;
end
Lat = -90;
for i = 2:NLat,
    Lat(i) = Lat(i-1) + 0.5;
end

%%
GDL = shaperead('GDL Shapefiles V4.shp');

%%
for i = 1:size(GDL),
    [z, r]= vec2mtx(GDL(i).Y,GDL(i).X,Res,[-90 90],[-180 180]); %z == 1 if border
    y = z; %fill in
    for j=2:359,
        for k=2:719,
            if z(j,k)==0 & sum(z(1:j-1,k))>0 & sum(z(j+1:360,k))>0 & sum(z(j,1:k-1))>0 & sum(z(j,k+1:720))>0
                y(j,k) = 1;
            end
        end
    end
    T0 = B.*y;
    T(i) = sum(nansum(T0))/sum(sum(y));
end