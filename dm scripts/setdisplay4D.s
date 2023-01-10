image f:=getfrontimage()
imagedisplay fd= f.imagegetimagedisplay(0)
number x,y
		x = f.imagegetdimensionsize(2)
		y = f.imagegetdimensionsize(3)
			
imagedisplaysetdisplayedlayers(fd, 0, x-1, 0, y-1)