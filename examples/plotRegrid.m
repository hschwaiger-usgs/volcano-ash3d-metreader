clear all;

dat1 = load('out02_1_2.dat'); % CompH
%dat1 = load('out05_2_2.dat'); % CompP
%dat1 = load('out05_3_2.dat'); % MetH
%dat1 = load('out05_4_2.dat'); % MetP

i1   = dat1(:,1);
j1   = dat1(:,2);
k1   = dat1(:,3);
x1   = dat1(:,4);
y1   = dat1(:,5);
z1   = dat1(:,6);
v1   = dat1(:,7);

nx0 = min(i1);
ny0 = min(j1);
nz0 = min(k1);

nx1 = max(i1);
ny1 = max(j1);
nz1 = max(k1);

if (nz0==nz1)
  % map view (z=something)
  [nx1 ny1 nz1]
  X1 = reshape(x1,ny1,nx1);
  Y1 = reshape(y1,ny1,nx1);
  V1 = reshape(v1,ny1,nx1);
  figure;
  colormap(jet);
  imagesc(X1,Y1,V1);axis xy
elseif (nx0==nx1)
  % Meridional transect (x=something)
  [nx1 ny1 nz1]
  Z1 = reshape(z1,nz1,ny1);
  Y1 = reshape(y1,nz1,ny1);
  V1 = reshape(v1,nz1,ny1);
  figure;
  colormap(jet);
  imagesc(Y1,Z1,V1);axis xy
elseif (ny0==ny1)
  %  Zonal transect (y=something)
  [nx1 ny1 nz1]
  Z1 = reshape(z1,nz1,nx1);
  X1 = reshape(x1,nz1,nx1);
  V1 = reshape(v1,nz1,nx1);
  figure;
  colormap(jet);
  imagesc(X1,Z1,V1);axis xy
  %contour(X1,Z1,V1)
end
