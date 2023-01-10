// $BACKGROUND$

image f:=getfrontimage()
imagedisplay df = f.imagegetimagedisplay(0)
number kx , ky,i,j, l,r,b,t
number x, y 


kx = f.ImageGetDimensionSize(0)
ky = f.ImageGetDimensionSize(1)
x = f.ImageGetDimensionSize(2)
y = f.ImageGetDimensionSize(3)

roi nroi = newroi()
nroi = df.imagedisplaygetroi(0)
nroi.ROIgetrectangle(t,l,b,r)

result("t "+t+" l "+l+" b "+b+" r "+r+" \n")
image out := realimage("",4,r-l,b-t, x,y)
//image out := realimage("",4,kx,ky)


	out = f.SliceN(4,4, l,t,0,0,  0,r-l,1,   1,b-t,1, 2,x,1, 3,y,1)

showimage(out)
imagedisplay do = out.imagegetimagedisplay(0)
do.imagedisplaysetdisplayedlayers(0,x-1,0,y-1)

