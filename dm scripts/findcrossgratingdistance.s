// $BACKGROUND$

// Script to calculate the average distance of the squares in a cross grating sample image.
// Just make sure, that the cross grating image is the front most image and execute the script

// It calculates the maximas of the autocorrelation image and the puts this pixel and its surroundings to zero
// and then searches for the next maxima.
// It assumes that the old calibration is somewhat correct. 



// subpixel info from V. Hou
Number ImageRefineExtrema(Image img, number &px, number &py)
{

number tx1, tx2, ty1, ty2
number numa, numb, numc
number valX, valY

// do a parabolic fit in direction x
numa = img.GetPixel(px-1, py)
numb = img.GetPixel(px , py)
numc = img.GetPixel(px+1, py)
tx2=(numa-numc)/2.0
tx1=(numa+numc)/2.0-numb
valX=(2*numa*numC-numa**2-numc**2)/8/(numa+numc-2*numb)+numb

// do a parabolic fit in direction y
numa = img.GetPixel(px, py-1)
numb = img.GetPixel(px, py )
numc = img.GetPixel(px, py+1)
ty2=(numa-numc)/2.0
ty1=(numa+numc)/2.0-numb
valY=(2*numa*numC-numa**2-numc**2)/8/(numa+numc-2*numb)+numb

// update pixel position
px += tx2/(2.0*tx1)
py += ty2/(2.0*ty1)

// check to see whether to return mimimum or maximum
If(numa>=numb && numc>=numb) return min(min(valX, valY), numb)
If(numa<=numb && numc<=numb) return max(max(valX, valY), numb)

Return numb
}


// get front image info
image n,nm,m,nn,his,l,ll,jj
nm := getfrontimage()
n = IFMApplyFilter( nm, "Hanning Window (default)")
number scale,org,x,y
string units
getsize(n,x,y)

ImageGetDimensionCalibration( nm, 0,org, scale, units,0)


// size of the cross grating distance in nm
number lent = 462.963
if(units != "nm") lent/=1000

// size of the cross grating distance in pixels
number pre = lent/scale

// half size of the square to be cut away after the maxima has been found
number cut= pre/2 +5

// estimated number of maximas in the image
number num= floor(x/pre)**2



number pycmi,pxcmi,pycma,pxcma,j


//IFMApplyFilterInPlace( m,"Combined Filter 0" )
nn := AutoCorrelate( n)

m=nn


number i,max,px,py,min,sum,mo1,mo2,c=2,d=2
number one,two,three,four,five,six,ppx,ppy,av
image cx = exprsize(num,1,0)
image cy = exprsize(num,1,0)
image tt = exprsize(num,1,0)
image ave = exprsize(num,1,0)

number num1=0


for(i=0;i<num;i++)
{
max = max(m,px,py)
openandsetprogresswindow(""+(((i+1)/num)*100)+"%","phase 1","")

	if(px>1 && px<x-2 && py>1 && py< y-2)
	{
	
	// find maxima and save the pixel coordinates
	ImageRefineExtrema(m, px, py)

	cx[i,0] = px
	cy[i,0] = py

	num1++
	}
	
	// make sure that the cut does not go over image boundaries
	pycmi = py-cut
	pxcmi = px-cut
	pycma = py+cut
	pxcma = px+cut
    
	if (pycmi<0) pycmi=0
	if (pxcmi<0) pxcmi=0
	if (pycma>y-1) pycma=y
	if (pxcma>x-1) pxcma=x

    // cut the maximas surroundings
	m[pycmi,pxcmi,pycma,pxcma] = 0
	
//showimage(m)
//sleep(0.5)
}

// get the new number of maximas
num = num1

// image for the distances
image k = exprsize(num,num,0)


for(i=0;i<num;i++)
{
px = getpixel(cx,i,0)
py = getpixel(cy,i,0)

// calculate the distances of every pixel to every pixel, here the diagonal will be zero (points distance to itself)
// and every distance will be counted twice (mirror symmetry across the diagonal). 
// Not perhaps the most fastest way to do this...
k[icol,i] += sqrt((px - cx)**2 + (py - cy)**2) 

}

number sum1,aa,count=0

// choose the values around the pre value and put them to value 1. This depends on the old calibration, if something is going wrong, it is here
// this is for the nubmer of samples
ll = tert(k> (pre*(8/10)) && k<(pre*(12/10)),1,0)


// sum is the number of samples (distances) calculated. (real one is this divided by two)
sum = sum(ll)


// choose the values around the pre value. This depends on the old calibration, if something is going wrong, it is here
// this is for the nubmer of samples
ll = tert(k> (pre*(8/10)) && k<(pre*(12/10)),k,0)

// the average distance
sum1 = sum(ll)/sum


// Do a check image
image check = exprsize(sum,1,0)
image check1 = exprsize(sum,1,0)
image check2 = exprsize(sum,1,0)

// go through the whole image
for(i=0;i<num;i++)
{
	for(j=0;j<num;j++)
	{
	aa = getpixel(ll,i,j)
		if(aa!=0) 
		{
        
        // the distance - average distance
        check1[count,0] = aa 
		check[count,0] = aa -sum1
		count++
		}

	}
openandsetprogresswindow(""+(((i+1)/num)*100)+"%","phase 2","")
}





showimage(check)
setname(check, "Here 0 is the mean size, is there too many outliers???")
lineplotimagedisplay ldisp = imagegetimagedisplay(check,0)
LinePlotImageDisplaysetSliceDrawingStyle( ldisp,0,1)

number maxc,res,count1=0,summ
maxc = max(check)
GetNumber( "Give a abs(number) to cut away outliers, or if all is ok, just press ok ",maxc+10, res )

// count the outliers
check2 = tert(abs(check) > res, 1,0)
summ = sum(check2)

image check3 = exprsize(sum-summ,1,0)


// create the final image without outliers
for(i=0;i<sum;i++)
{
	if(getpixel(check2,i,0)==0)
	{
	check3[count1,0] = getpixel(check1,i,0)
	count1++
	}
}


deleteimage(check)


sum1 = sum(check3)/(sum-summ)
//check3-=sum1
if(units!="nm") check3 = 1000*check3*(lent/sum1)/sum1
else  check3 = check3*(lent/sum1)/sum1

showimage(check3 )
setname(check3, "final pixel sizes")
lineplotimagedisplay ldisp1 = imagegetimagedisplay(check3,0)
LinePlotImageDisplaysetSliceDrawingStyle( ldisp1,0,1)

av = sqrt(variance(check3))

result("\n\n")
result("average length of the squares in pixels is "+sum1+" \n") 
result("In old calibrations the length is "+(scale*sum1)+" "+units+" \n")
result("Old Pixel size is "+scale+" "+units+"\n")
result("Pixel size according to this calculation is "+(lent/sum1)+" "+units+" with st_dev "+av+" calculated from "+((count1/2))+" samples (distances)\n")


if(OkCancelDialog( "Do you want to chance the calibration of the image to the new value?" )==1) setscale(n,(lent/sum1),(lent/sum1))




