image n:=getfrontimage()
imagedisplay nd = n.imagegetimagedisplay(0)
roi mroi = nd.imagedisplaygetroi(0)
number t,l,b,r
number kx = n.imagegetdimensionsize(2)
number ky = n.imagegetdimensionsize(3)
number sx = n.imagegetdimensionscale(0)
number sy = n.imagegetdimensionscale(1)
mroi.roigetrectangle(t,l,b,r)
string name,sc
sc =  n.ImageGetDimensionUnitString( 0 )
getname(n,name)
taggroup ntag = n.imagegettaggroup()
image new := realimage("cut_" + name, 4,r-l,b-t, kx,ky)
taggroup newtag = new.imagegettaggroup()
new = n.slicen(4,4,  l,t,0,0,  0, r-l,1,  1,b-t,1,  2,kx,1,  3,ky,1)
new.imagesetdimensionscale(0,sx)
new.imagesetdimensionscale(1,sy)
taggroupcopytagsfrom(newtag,ntag)
new.ImagesetDimensionUnitString( 0,sc )
new.ImagesetDimensionUnitString( 1,sc )

showimage(new)

imagedisplay dnew = new.imagegetimagedisplay(0)

new.setsurveytechnique(2)
dnew.ImageDisplaySetDoAutoSurvey(1 )
dnew.ImageDisplaySetDisplayedLayers( 0,kx-1,0,ky-1 )
