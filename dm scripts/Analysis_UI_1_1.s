// $BACKGROUND$

//image fro  := OpenImage(  "C:/Jup/Hasan_orig/4D STEM2.dm4" )
//showimage(fro)
number roion=0,threshold=-1


// The Class for the dialog

Class Dialog_UI : UIFrame
	{
	
	Number ImageRefineExtrema(object self,Image img, number &px, number &py)
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
	image SU_GetCC(object self, number &px, number &py, image src, image mov,  number &qval, number mode, number refine) //microns
		{

		number mi,dimension, origin, scale,calibrationFormat

		number width, height
		string unit
		getsize(src,width, height)
		ImageGetDimensionCalibration( src, dimension, origin, scale,unit, calibrationFormat )					
		Image xcorr1 = CreateFloatImage("Xc",width, height);
			
			if(mode==0 )
			{
			image prevF = IFMapplyfilter(mov, "AutoEM")
			image srcF = IFMapplyfilter(src, "AutoEM")
			xcorr1 = srcF.CrossCorrelate(prevF)
			}
			if(mode==1)
			{
			xcorr1 = src.CrossCorrelate(mov)
			}
		
			// find the maximum and its coordinates in correlation image
		Number spotX, spotY
		qval = max(xcorr1, spotX, spotY)
		
		if(refine) self.ImageRefineExtrema(xcorr1 ,spotx, spoty)
		px = ( width/2 - spotx) ;  py = ( height/2 - spoty)
		//result("px "+px+" py "+py+"  \n");
		//showimage(xcorr1)
		return xcorr1
		}

	void findcentreofgravity(object self,image centralroi, number &cogx, number &cogy)
		{
		// For information on this calculation see John Russ - Image Processing Handbook, 2nd Ed. p489

			number xpos, ypos, maxval, imgsum, xsize, ysize, i
			image tempimg

			maxval=max(centralroi, xpos, ypos)
			imgsum=sum(centralroi)
	
			getsize(centralroi, xsize, ysize)

			// Traps for a blank image
			
			if(imgsum==0) // the image is blank so set the CoGs to the centre of the image and return
				{
					cogx=(xsize-1)/2 //minus one since the centre of a 2 x 2 image is 0.5,0.5 - 0,0 is a position
					cogy=(ysize-1)/2 //minus one since the centre of a 2 x 2 image is 0.5,0.5 - 0,0 is a position
					return
				}


			// Collapse the image down onto the x axis

			image xproj=realimage("",4,xsize,1)
			xproj[icol,0]+=centralroi


			// Rotate the passed in image through 90 degs so that rows become columns
			// Then collapse that image down onto the x axis (was the y axis)
			
			tempimg=realimage("", 4, ysize, xsize)
			tempimg=centralroi[irow,icol]
			image yproj=realimage("",4,ysize, 1)
			yproj[icol,0]+=tempimg 

			yproj=yproj*(icol+1) // NB the +1 ensures that for the left column and top row, where
			xproj=xproj*(icol+1) // icol=0 are included in the weighting. 1 must be subtracted from
								// the final position to compensate for this shift
			cogx=sum(xproj)
			cogy=sum(yproj)
			cogx=(cogx/imgsum)-1 - xsize/2 // compensation for the above +1 to deal with icol=0
			cogy=(cogy/imgsum)-1 - ysize/2
			if(mod(xsize,2)==0) cogx = cogx + 0.5
			if(mod(ysize,2)==0) cogy = cogy + 0.5
		}
	void get_CoM_images(object self,image datacube, image mask, number normalize, image &CoMx, image &CoMy)
		{
		number x,y,i,j,kx,ky
		
		x = datacube.imagegetdimensionsize(0)
		y = datacube.imagegetdimensionsize(1)
		kx = datacube.imagegetdimensionsize(2)
		ky = datacube.imagegetdimensionsize(3)
		
		image DP:=realimage("", 4, kx,ky)
		image mass:=realimage("", 4, kx,ky)
		
			for(j=0;j<y;j++)
			{
				for(i=0;i<x;i++)
				{
				DP = slicen(datacube, 4,2, i,j,0,0, 2,kx,1, 3,ky,1) 
				mass[j,i] = sum(DP*mask)
				CoMx[j,i] = sum((irow+1)*DP) / mass[j,i]
				CoMy[j,i] = sum((icol+1)*DP) / mass[j,i]
				}
			}
			
			if(normalize)
			{
			CoMx -= mean(CoMx)
			CoMy -= mean(CoMy)
			}
		}
		
	void AddROIbuttonresponse(object self )
		{
		TagGroup PTag = GetPersistentTagGroup( )
		number x,y,nr=0,addROIC,addROIf,i, cubeID, planeID
	
		image front := getfrontimage()
		taggroup fronttag = front.imagegettaggroup()
		cubeID = front.imagegetid()
		
		x = front.imagegetdimensionsize(0)
		y = front.imagegetdimensionsize(1)
		
		imagedisplay idisp = imagegetimagedisplay(front,0)
		nr = idisp.imagedisplaycountrois()

		if(nr==0) TagGroupSetTagAsBoolean(PTag, "SU:ROIListenerOn", 0);
		self.dlggetvalue("addROICheckBox", addROIC)
		self.dlggetvalue("addROIfield", addROIf)
		if(addROIf<1) addROIf=1; if(addROIf>63) addROIf=63;
		if(nr>0)
		{
		roi oroi
		for(i=0;i<nr ;i++)
		{
		oroi = idisp.imagedisplaygetroi(i)
		if(ImageDisplayIsROISelected( idisp, oroi )) break
		}
		
		string lab = oroi.roigetlabel()
			for(i=0;i<addROIf ;i++)
			{
			roi nroi = newroi()	
			nroi.roisetrectangle( round(y/2)-1, round(x/2)-(addROif - i ),round(y/2),round(x/2))
						
			if(addROIC) ROISetLabel( nroi, ""+(nr+1+i)+"" )
			else ROISetLabel( nroi, lab )
			idisp.imagedisplayaddroi(nroi)
			}
		}
		else
		{
			
			for(i=0;i<addROIf ;i++)
			{
			roi nroi = newroi()	
			nroi.roisetrectangle( round(y/2)-1, round(x/2)-(addROif - i ),round(y/2),round(x/2))
						
				if(addROIC) ROISetLabel( nroi, ""+(nr+1+i)+"" )
				else 
				{
				if(nr==0)nr=1
				ROISetLabel( nroi, ""+nr+"" )	
				}
				idisp.imagedisplayaddroi(nroi)
			}
		}
		TagGroupSetTagAsBoolean(PTag, "SU:defaults:addRoiCheckBox", addROIC);
		TagGroupSetTagAsUInt32(PTag,  "SU:defaults:addRois", addROIf);
		TagGroupGetTagAsBoolean(PTag, "SU:ROIListenerOn", roion);
		sleep(0.2)
		if(!roion) {SU_ROIThread(front);}

		return;
		}
	void AddROIbuttonresponse1(object self )
		{
		TagGroup PTag = GetPersistentTagGroup( )
		number x,y,nr=0,max,px,py,addROIC,addROIf,i
	
		image front := getfrontimage()
		max = max(front, px,py)

		x = front.imagegetdimensionsize(0)
		y = front.imagegetdimensionsize(1)
		
		imagedisplay idisp = imagegetimagedisplay(front,0)
		nr = idisp.imagedisplaycountrois()
		
		self.dlggetvalue("addROICheckBox", addROIC)
		self.dlggetvalue("addROIfield", addROIf)
		
		if(addROIf<1) addROIf=1; if(addROIf>63) addROIf=63;
		
		if(nr>0)
		{
		roi oroi
		for(i=0;i<nr ;i++)
		{
		oroi = idisp.imagedisplaygetroi(i)
		if(ImageDisplayIsROISelected( idisp, oroi )) break
		}
		string lab = oroi.roigetlabel()
		
			for(i=0;i<addROIf ;i++)
			{
			roi nroi = newroi()	
			nroi.roisetoval( py-9, px-9 -(addROif - i ) ,py+9,px+9 -(addROif - i ))
			if(addROIC) ROISetLabel( nroi, ""+(nr+1+i)+"" )
			else ROISetLabel( nroi, lab )
			idisp.imagedisplayaddroi(nroi)
			}
		}
		else
		{
			
			for(i=0;i<addROIf ;i++)
			{
			roi nroi = newroi()	
			nroi.roisetoval( py-9, px-9 -(addROif - i ) ,py+9,px+9 -(addROif - i ))
						
				if(addROIC) ROISetLabel( nroi, ""+(nr+1+i)+"" )
				else 
				{
				if(nr==0){ nr=1;}
				ROISetLabel( nroi, ""+nr+"" )	
				}
			idisp.imagedisplayaddroi(nroi)
			}
		}
		TagGroupSetTagAsBoolean(PTag, "SU:defaults:addRoiCheckBox", addROIC);
		TagGroupSetTagAsUInt32(PTag,  "SU:defaults:addRois", addROIf);
		TagGroupGetTagAsBoolean(PTag, "SU:ROIListenerOn", roion);
		sleep(0.2)
		if(!roion) {SU_ROIThread(front);}
	
		return;
		}
			
	void addVHAADFbuttonResponse( object self)
		{
		TagGroup PTag = GetPersistentTagGroup( )
		number VHAADFin, VHAADFout, NSegments, NSegmentsrot,ison,s, cubeid, planeid
		
		self.dlggetvalue("VHAADFinfield",VHAADFin)
		self.dlggetvalue("VHAADFoutfield", VHAADFout)
		self.dlggetvalue("VHAADFSegsfield",NSegments)
		self.dlggetvalue("VHAADFSegsRotfield",NSegmentsrot)
	
		TagGroupSetTagAsDouble(PTag, "SU:defaults:VHAADFin",VHAADFin);
		TagGroupSetTagAsDouble(PTag, "SU:defaults:VHAADFout",VHAADFout);
		TagGroupSetTagAsuint16(PTag, "SU:defaults:NSegments", NSegments);
		TagGroupSetTagAsDouble(PTag, "SU:defaults:NSegmentsrot",NSegmentsrot);
	
		number x,y,nr=0,lx1,lx2,ly1,ly2,i,j,mx,my,maxIm,delta=0, deltastep=0, dlx1,dlx2,dly1,dly2,nim,kx,ky
		string type
		image cube, plane
		nim = countimages()
		image dummy := getfrontimage()
	

		for(i=0;i<nim;i++)
		{
			if(dummy.imageisvalid())
			{

			taggroup ftags = dummy.imagegettaggroup()
			TagGroupGetTagAsString( ftags,"SU:Auxiliary:type", type)
			if(type == "cube") {cubeid = imagegetid(dummy);}
			if(type == "plane") {planeid = imagegetid(dummy);}//result("planeid "+planeid+" \n");
			type=""
			dummy := findnextimage(dummy)
			}
	
		}
		cube := findimagebyid(cubeid)
		kx = cube.imagegetdimensionsize(2)
		ky = cube.imagegetdimensionsize(3)
		plane := findimagebyid(planeid)
		if(VHAADFin>(1+sqrt(2)*kx/2) || VHAADFout>(1+sqrt(2)*kx/2) ) {result("ROI is bigger than image (1+ sqrt(2)*x_halfsize): "+(1+sqrt(2)*kx/2)+" \n"); return;}
		if(VHAADFin>(1+sqrt(2)*ky/2) || VHAADFout>(1+sqrt(2)*ky/2) ) {result("ROI is bigger than image (1+ sqrt(2)*y_halfsize): "+(1+sqrt(2)*ky/2)+" \n"); return;}
		
		s = ImageGetNumDimensions(cube)	
		x = plane.imagegetdimensionsize(0)
		y = plane.imagegetdimensionsize(1)
		
		if(kx!=x){result("Image "+getname(plane)+" , id "+imagegetid(plane)+" size "+x+" does not match with 4D image "+kx+" \n"); return;}

		//maxIm=max(plane,mx,my)
		mx = x/2 ; my = y/2;
		//maxIm=max(plane,mx,my)
		if(mod(x,2)==0) mx = mx - 0.5
		if(mod(y,2)==0) my = my - 0.5
		imagedisplay idisp = imagegetimagedisplay(plane,0)
		nr = idisp.imagedisplaycountrois()
		
		if(nr>0) for(i=0;i<nr;i++) idisp.imagedisplaydeleteroi(ImageDisplayGetROI(idisp, 0 ))
		if(nr==0) TagGroupSetTagAsBoolean(PTag, "SU:ROIListenerDFOn", 0);
			
		nr = plane.CountAnnotations()
		if(nr>0) for(i=0;i<nr;i++) plane.DeleteAnnotation( GetNthAnnotationID( plane,0 ))

		roi nroi_in = newroi()
		if(VHAADFin>0)	
		{
		nroi_in.roisetoval( my-VHAADFin, mx-VHAADFin, my+VHAADFin, mx+VHAADFin)
		TagGroupSetTagAsLong(PTag,  "SU:Auxiliary:VHAADFin_ID", ROIGetID(nroi_in ));
		idisp.imagedisplayaddroi(nroi_in)
		}
		else TagGroupSetTagAsLong(PTag,  "SU:Auxiliary:VHAADFin_ID", 0);
	
		roi nroi_out = newroi()	
		nroi_out.roisetoval( my-VHAADFout, mx-VHAADFout, my+VHAADFout, mx+VHAADFout)
		TagGroupSetTagAsLong(PTag,  "SU:Auxiliary:VHAADFout_ID", ROIGetID(nroi_out ));
		idisp.imagedisplayaddroi(nroi_out)
		
			if(NSegments>0)
			{
			deltastep = 180/ (NSegments)
				for(i=0;i< NSegments;i++)
				{
				delta = NSegmentsrot*(pi()/180)  + i*deltastep*(pi()/180)
				
				roi lineroi = newroi()
				lx1 = 0
				ly1 = -VHAADFout
				lx2 = 0
				ly2 = VHAADFout
				dlx1 = mx +cos(delta)*lx1 - sin(delta)*ly1;  dly1 = my + sin(delta)*lx1 + cos(delta)*ly1;
				dlx2 = mx + cos(delta)*lx2 - sin(delta)*ly2;  dly2 = my + sin(delta)*lx2 + cos(delta)*ly2;

				roiSetLine(lineroi,dlx1,dly1,dlx2,dly2)
				TagGroupSetTagAsUInt32(PTag,  "SU:Auxiliary:VHAADFline_ID_"+(i+1)+"", ROIGetID(lineroi ));
				TagGroupSetTagAsdouble(PTag,  "SU:Auxiliary:VHAADFline_ID_"+(i+1)+":alpha", (atan( (dly2 -my) / (dlx2 -mx) ) ));
				idisp.imagedisplayaddroi(lineroi)

				}
			deltastep = 360/ (NSegments*2)	
				for(i=0;i< NSegments;i++)
				{
				number rad =  VHAADFin + (VHAADFout  )/2
					for(j=0;j< NSegments*2;j++)
					{
					
					delta = NSegmentsrot*(pi()/180)  + j*deltastep*(pi()/180)
					
					if(j*deltastep>179 && NSegments>1) delta = delta +20*(pi()/180)
					if(NSegments==1){lx1 = rad; ly1 = 0; }
					else{lx1 = 0; ly1 = -rad; }

					dlx1 = mx +cos(delta)*lx1 - sin(delta)*ly1;  dly1 = my + sin(delta)*lx1 + cos(delta)*ly1;

					component tcomp = NewTextAnnotation( dlx1,dly1,""+(j+1)+"", x/20)
					idisp.componentaddchildatend(tcomp ) 
					tcomp .componentsetforegroundcolor(1,0,1) // sets the foreground colour to magenta
					}
				}
			}
		
		TagGroupGetTagAsBoolean(PTag, "SU:ROIListenerDFOn", ison);
		if(!ison) SU_ShowVirtualDarkFieldImage(cube,plane)
		}
	void alignbuttonResponse( object self)
		{	
		number val, diameter, boxsize
		self.dlggetvalue("alignpopup",val)
		self.dlggetvalue("alignfield",diameter)
		
		if(diameter<6000) boxsize = 8192
		if(diameter<3200) boxsize = 4096
		if(diameter<1600) boxsize = 2048
		if(diameter<800) boxsize = 1024
		if(diameter<400) boxsize = 512
		if(diameter<200) boxsize = 256
		if(diameter<70) boxsize = 128
		if(diameter<30) boxsize = 64
		if(diameter<10) boxsize = 32
		if(diameter<5) boxsize = 16
		
		number x,y,kx,ky,i,j,px,py,q, mx,my,mmx,mmy
		image f := getfrontimage()
		image nm,nn
		x = f.imagegetdimensionsize(0)
		y = f.imagegetdimensionsize(1)
		kx = f.imagegetdimensionsize(2)
		ky = f.imagegetdimensionsize(3)

		for(i=0;i<x;i++)
		{
		OpenAndSetProgressWindow( "St4DeM", "Aligning..", " "+ (i+1) +" / "+x+" ")
			for(j=0;j<y;j++)
			{
			nm = slicen(f, 4,2, i,j,0,0, 2,kx,1, 3,ky,1) 
			if(val==1)
			{
			nn := AutoCorrelate( nm)**2
			self.SU_GetCC(px, py, nm, nn,q,0,1) 
			}
			else if(val==2)self.findcentreofgravity(nm,px,py)
			else if(val==3)
			{
			q = max(nm, mmx,mmy)
			Image parsToFit := NewImage("tmp", 2, 6)
			parsToFit = 1
			Number chiSqr = 1e6
			Number conv_cond = 0.00001

			image Imbox = nm[mmy- boxsize/2, mmx- boxsize/2, mmy + boxsize/2, mmx + boxsize/2]
			Image errors := Imbox.ImageClone()
			errors = tert(abs(Imbox) > 1, sqrt(abs(Imbox)), 1)
			q = max(Imbox, mx,my)
			
			Image pars := NewImage("Gaussian2D Pars", 2, 6)
			pars = diameter
			pars[0,0] = q  // estimate normalization from peak of data
			pars[5,0] = 0   // 100 radians doesn't make sense
			pars[1,0] = mx     // center in x
			pars[3,0] = my     // center in y
			
			Number ok = FitGaussian2D(Imbox, errors, pars, parsToFit, chiSqr, conv_cond)
			//px = pars[1,0]
			//py = pars[3,0]
			px = kx/2 -( mmx - boxsize/2 + getpixel(pars,1,0) ) ; 
			py = ky/2 -( mmy - boxsize/2 + getpixel(pars,3,0) )
			}
			else return
			//max(nm,px,py)
			//px = px - kx/2
			//py = py - ky/2
			slicen(f, 4,2, i,j,0,0, 2,kx,1, 3,ky,1) = warp(nm, icol + px, irow + py)
			
			}
		//OpenAndSetProgressWindow( ""+i+"/"+y+"", "", "")
		//result("x: "+i+" \n")
		}
		OpenAndSetProgressWindow( "", "", "")

		}
		
	void alignbuttonResponseT( object self)
		{
		self.startthread("alignbuttonResponse");
		}
	void VHAADFImShowbutton( object self )
		{
		TagGroup PTag = GetPersistentTagGroup( )
		number NSegments,alpha= - pi(), inid,outid,doin=0,VHAADFin,VHAADFout,NSegmentsrot,maxIm,mx,my,x,y,i
		number cubeid, planeid,VHAADFIm1,VHAADFIm2,nim,a1,a2,a3,a4,kx,s
		roi roi_in,roi_out
		
		self.dlggetvalue("VHAADFinfield",VHAADFin)
		self.dlggetvalue("VHAADFoutfield", VHAADFout)
		self.dlggetvalue("VHAADFSegsfield",NSegments)
		self.dlggetvalue("VHAADFSegsRotfield",NSegmentsrot)
		self.dlggetvalue("VHAADFIm1field",VHAADFIm1)
		self.dlggetvalue("VHAADF2m1field",VHAADFIm2)
			
		NSegmentsrot = (NSegmentsrot*pi())/180
		number anglestep = pi()/NSegments
		
		TagGroupGetTagAsLong(PTag,  "SU:Auxiliary:VHAADFin_ID", inid);
		if(inid>0) {doin=1;  roi_in = GetROIFromID( inid );}
		TagGroupGetTagAsLong(PTag,  "SU:Auxiliary:VHAADFout_ID", outid);
		roi_out = GetROIFromID( outid );

		string type
		image cube, plane
		nim = countimages()
		image dummy
		dummy.getfrontimage()

			for(i=0;i<nim;i++)
			{
				if(dummy.imageisvalid())
				{

				taggroup ftags = dummy.imagegettaggroup()
				TagGroupGetTagAsString( ftags,"SU:Auxiliary:type", type)
				if(type == "cube") {cubeid = imagegetid(dummy);}
				if(type == "plane") {planeid = imagegetid(dummy);}
				type=""
				dummy := findnextimage(dummy)
				}
			
			}

		cube := findimagebyid(cubeid)
		kx = cube.imagegetdimensionsize(2)
		image mask = findimagebyid(planeid)
		
		s = ImageGetNumDimensions(cube)	
		x = mask.imagegetdimensionsize(0)
		y = mask.imagegetdimensionsize(1)

		if(kx!=x){result("Image "+getname(plane)+" , id "+imagegetid(dummy)+" size "+x+" does not match with 4D image "+kx+" \n"); return;}
		
		maxIm=max(mask,mx,my)
		x = mask.imagegetdimensionsize(0)
		y = mask.imagegetdimensionsize(1)
		
		mask = 0
		image mask1= mask
		
		AddArrayMaskToImage(mask, x/2, y/2, VHAADFout, VHAADFout, 0, 0, x, y,0,0 )
		if(doin){ AddArrayMaskToImage(mask1, x/2, y/2, VHAADFin, VHAADFin, 0, 0, x, y,0,0 ); mask = mask - mask1; mask1 = mask;}

			for(i=0;i< NSegments*2;i++)
			{
			mask = tert( ( (mask ==1) && ( itheta >= alpha  ) && ( itheta < alpha + anglestep )   ), i+1, mask) //&& ( itheta < 0  )
			alpha = alpha + anglestep
			}
		
		if(!NSegmentsrot){RotateRight( mask); }//showimage(mask);
		else 
		{
		mask1 := rotate(mask, -pi()/2 - NSegmentsrot);
		number x1,y1
		getsize(mask1, x1,y1)
		mask := mask1[y1/2 -y/2, x1/2 - x/2, y1/2+y/2, x1/2 + x/2 ]
		//showimage(mask);
		}
		
		SU_ShowDPCImageThread(cube,mask,VHAADFIm1,VHAADFIm2,NSegmentsrot)
		}
	
	void iDPCImShowbuttonResponse( object self )
		{
		image x,y
		number kx,ky
		gettwoimages("x,y",x,y)
		getsize(x,kx,ky)

		compleximage aDPC = compleximage("aDPC", 8, kx,ky)
		setname(aDPC, "aDPC")
		real(aDPC) = x
		imaginary(aDPC) = y
		showimage(aDPC)

		image iDPC

		iDPC = realifft(realfft(x) + realfft(y))
		setname(iDPC, "iDPC")
		showimage(iDPC)

		}
		
	TagGroup CreateDialog_UI( object self )
		{
		TagGroup PTag = GetPersistentTagGroup( )
		number VHAADFin=10.0,VHAADFout=255.0,NSegments=1,NSegmentsrot=0.0, VHAADFIm1=1, VHAADFIm2=2,addRoiCheck=0, AddRois=1
		string version = "null"

		if (TagGroupDoesTagExist(PTag,"SU:defaults:VHAADFin") == 0) 	TagGroupSetTagAsDouble(PTag, "SU:defaults:VHAADFin",VHAADFin);
		if (TagGroupDoesTagExist(PTag,"SU:defaults:VHAADFout") == 0) 	TagGroupSetTagAsDouble(PTag, "SU:defaults:VHAADFout",VHAADFout);
		if (TagGroupDoesTagExist(PTag,"SU:defaults:NSegments") == 0) 	TagGroupSetTagAsuint16(PTag, "SU:defaults:NSegments", NSegments);
		if (TagGroupDoesTagExist(PTag,"SU:defaults:NSegmentsrot") == 0) 	TagGroupSetTagAsDouble(PTag, "SU:defaults:NSegmentsrot",NSegmentsrot);
		if (TagGroupDoesTagExist(PTag,"SU:defaults:VHAADFIm1") == 0) 	TagGroupSetTagAsuint16(PTag, "SU:defaults:VHAADFIm1", VHAADFIm1);
		if (TagGroupDoesTagExist(PTag,"SU:defaults:VHAADFIm2") == 0) 	TagGroupSetTagAsuint16(PTag, "SU:defaults:VHAADFIm2", VHAADFIm2);
		

		if (TagGroupDoesTagExist(PTag,"Private:Configuration:ApplicationVersion_2")) {TagGroupGetTagAsString(PTag,"Private:Configuration:ApplicationVersion_2", version);}
		if(version[0] == 49) TagGroupSetTagAsUInt32(PTag,"SU:DMVersion", 1);
		if(version[0] == 50) TagGroupSetTagAsUInt32(PTag,"SU:DMVersion", 2);
		if(version[0] == 51) TagGroupSetTagAsUInt32(PTag,"SU:DMVersion", 3);
		TagGroupSetTagAsBoolean(PTag, "SU:ROIListenerOn", 0);

		TagGroupGetTagAsDouble(PTag, "SU:defaults:VHAADFin",VHAADFin);
		TagGroupGetTagAsDouble(PTag, "SU:defaults:VHAADFout",VHAADFout);
		TagGroupGetTagAsuint16(PTag, "SU:defaults:NSegments", NSegments);
		TagGroupGetTagAsDouble(PTag, "SU:defaults:NSegmentsrot",NSegmentsrot);
		TagGroupGetTagAsuint16(PTag, "SU:defaults:VHAADFIm1", VHAADFIm1);
		TagGroupGetTagAsuint16(PTag, "SU:defaults:VHAADFIm2", VHAADFIm2);
		
		TagGroup Dialog_UI = DLGCreateDialog("4D-STEM")

				TagGroup tabs = DLGCreateTabList( 1 );

				TagGroup tabAnalysis = DLGCreateTab( "Analysis" ) 
				tabs.DLGAddTab(  tabAnalysis)		


				//Analysis
				taggroup addROIbutton = dlgcreatepushbutton("Add ROI: []","AddROIbuttonresponse").dlgidentifier("addROIbutton")
				taggroup addROIbutton1 = dlgcreatepushbutton("Add ROI: O","AddROIbuttonresponse1").dlgidentifier("addROIbutton1")
				taggroup addROICheckBox = DLGCreateCheckBox( "New Label",addRoiCheck).DLGIdentifier("addROICheckBox")
				taggroup addROIfield = DLGCreateintegerField(AddRois, 3).DLGIdentifier("addROIfield")
				taggroup analysisgroup = dlggroupitems(addROIfield, addROICheckBox, addROIbutton,addROIbutton1).dlgtablelayout(4,1,0)
				
				taggroup VHAADFinlabel = DLGCreateLabel("Radius Inner: ")	
				taggroup VHAADFoutlabel = DLGCreateLabel("Outer: ")
	
				taggroup VHAADF_box_items
				taggroup VHAADF_box = dlgcreatebox("Virtual DF/DPC", VHAADF_box_items).dlgexternalpadding(2,22).dlginternalpadding(12,12)
				taggroup VHAADFinfield = DLGCreateRealField( VHAADFin, 6,1).DLGIdentifier("VHAADFinfield")
				taggroup VHAADFoutfield = DLGCreateRealField( VHAADFout, 6,1).DLGIdentifier("VHAADFoutfield")
				taggroup VHAADFSegsfield = DLGCreateintegerField(NSegments, 3).DLGIdentifier("VHAADFSegsfield")
				taggroup VHAADFSegsRotlabel = DLGCreateLabel("Rotation:")
				taggroup VHAADFSegslabel = DLGCreateLabel("Lines: ")
				taggroup VHAADFSegsRotfield = DLGCreateRealField( NSegmentsrot, 5,1).DLGIdentifier("VHAADFSegsRotfield")
				taggroup VHAADFgroup1 = dlggroupitems(VHAADFinlabel, VHAADFinfield, VHAADFoutlabel,VHAADFoutfield).dlgtablelayout(4,1,0)
				taggroup VHAADFgroup2 = dlggroupitems(VHAADFSegslabel, VHAADFSegsfield,VHAADFSegsRotlabel,VHAADFSegsRotfield).dlgtablelayout(4,1,0)
				VHAADF_box_items.dlgaddelement(VHAADFgroup1)
				VHAADF_box_items.dlgaddelement(VHAADFgroup2)
				taggroup addVHAADFbutton = dlgcreatepushbutton("Add Mask: ","addVHAADFbuttonResponse")
				VHAADF_box_items.dlgaddelement(addVHAADFbutton)
				taggroup VHAADFImlabel = DLGCreateLabel(" - ")
				taggroup VHAADFIm1field = DLGCreateintegerField(VHAADFIm1, 2).DLGIdentifier("VHAADFIm1field")
				taggroup VHAADFIm2field = DLGCreateintegerField(VHAADFIm2, 2).DLGIdentifier("VHAADF2m1field")
				taggroup VHAADFImShowbutton = dlgcreatepushbutton("Show: ","VHAADFImShowbutton")
				taggroup VHAADFgroup3 = dlggroupitems(VHAADFIm1field,VHAADFImlabel, VHAADFIm2field,VHAADFImShowbutton ).dlgtablelayout(4,1,0)
				taggroup iDPCImShowbutton = dlgcreatepushbutton("Get iDPC","iDPCImShowbuttonResponse")
				
				VHAADF_box_items.dlgaddelement(VHAADFgroup3)
				VHAADF_box_items.dlgaddelement(iDPCImShowbutton)
				taggroup align_box_items
				taggroup align_box = dlgcreatebox("Align 4D",align_box_items).dlginternalpadding(12,12)
				taggroup alignpopup = DLGCreatePopup(2).DLGIdentifier("alignpopup")
				DLGAddPopupItemEntry( alignpopup, "AutoCor" ,0)
				DLGAddPopupItemEntry( alignpopup, "CoM" ,1)
				DLGAddPopupItemEntry( alignpopup, "GaussFit" ,2)
				
				taggroup alignbutton = dlgcreatepushbutton("Align","alignbuttonResponseT")
				taggroup alignfield = DLGCreateintegerField(5, 5).DLGIdentifier("alignfield")
				taggroup alignroup = dlggroupitems(alignbutton,alignpopup , alignfield).dlgtablelayout(3,1,0)
				align_box_items.dlgaddelement(alignroup)
				
				DLGAddElement(  tabAnalysis,  align_box)
				DLGAddElement(  tabAnalysis,  analysisgroup)
				DLGAddElement(  tabAnalysis,  VHAADF_box)
				 
				Dialog_UI.dlgaddelement(tabs)
				
				taggroup footer=dlgcreatelabel("T. Uusimaeki v1.2 Apr 2021").dlgexternalpadding(0,2)
				Dialog_UI.dlgaddelement(footer)

		return Dialog_UI
		}
		
		// default object constructor
		
		Dialog_UI( object self )
			{
		
				number dialogID=self.ScriptObjectGetID()
				//ScanIsOn = 0
				self.super.init( self.CreateDialog_UI() )
		
			}


		// default object destructor
		
		~Dialog_UI( object self )
			{
				//number dialogID=self.ScriptObjectGetID()
				//result("\nDialog with ID: "+dialogID+" destroyed.")
			}



	}

void main()
	{



		// Create the dialog
		object Dialog_UI = Alloc(Dialog_UI)
		
		Dialog_UI.Display("4D STEM")

		DocumentWindow window = Dialog_UI.GetFrameWindow();
		number appT, appL, appB, appR,winT, winL, winB, winR
		ApplicationGetBounds(appT, appL, appB, appR);
		window.WindowGetFrameBounds(winT, winL, winB, winR);
		window.WindowSetFramePosition(appR- (WinR-WinL), appT);

		Return
	}


main()




































