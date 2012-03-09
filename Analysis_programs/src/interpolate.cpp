#pragma once

#define debug 1

//return the three coefficient of the parabola passing for the three passed points
void parabolic_spline(double &a,double &b,double &c,double *xin,double *yd)
{
  double ya[3],yb[3],yc[3];
  for(int i=0;i<3;i++)
    {
      ya[i]=sqr(xin[i]);
      yb[i]=xin[i];
      yc[i]=1;
    }
  
  double D=det3(ya,yb,yc);
  double Da=det3(yd,yb,yc);
  double Db=det3(ya,yd,yc);
  double Dc=det3(ya,yb,yd);
  
  a=Da/D;
  b=Db/D;
  c=Dc/D;
}

bvec interpolate_multi(bvec xin,bvec yin,bvec xout,const char *outpath=NULL)
{
  bvec yout(xout.nel,xout.nboot,xout.njack);
  int nxin=yin.nel;
  int nxout=yout.nel;
  
  int nboot=yout.nboot;
  
  int nplot=100;
  double xplot[nplot];
  bvec yplot(nplot,xout.nboot,xout.njack);
  
  for(int ixout=0;ixout<nxout;ixout++)
    {
      //check if the needed x is not below the lower one or above the heavier one
      double dmin=(xin[1][nboot]-xin[0][nboot])/2;
      double dmax=(xin[nxin-1][nboot]-xin[nxin-2][nboot])/2;
      bool linbelow=(xout[ixout][nboot]-xin[0][nboot]<dmin);
      bool linabove=(xin[nxin-1][nboot]-xout[ixout][nboot]<dmax);
      if(nxin==2) linabove=1;
      
      for(int iboot=0;iboot<nboot+1;iboot++)
	//if linear
	if(linabove||linbelow)
	  {
	    //choose the pair
	    int sup;
	    if(linabove) sup=nxin-1;
	    else sup=1;
	    
	    //find the m and q
	    double m=(yin[sup][iboot]-yin[sup-1][iboot])/(xin[sup][iboot]-xin[sup-1][iboot]);
	    double q=yin[sup-1][iboot]-m*xin[sup-1][iboot];
	    
	    //interpolate
	    yout[ixout].data[iboot]=m*xout[ixout][iboot]+q;
	    
	    if(outpath!=NULL && ixout==0)
	      for(int iplot=0;iplot<nplot;iplot++)
		{
		  xplot[iplot]=xin[sup-1][iboot]+(xin[sup][iboot]-xin[sup-1][iboot])/99*iplot;
		  yplot[iplot].data[iboot]=xplot[iplot]*m+q;
		}
	  }
	else
	  {
	    //if parabolic interpolation, find the nearest point
	    int nearx=0;
	    for(int ixin=0;ixin<nxin;ixin++) if(fabs(xin[ixin][nboot]-xout[ixout][nboot])<fabs(xin[nearx][nboot]-xout[ixout][nboot])) nearx=ixin;
	    
	    //copy the x and y for current bootstrap
	    double txin[3]={xin[nearx-1][iboot],xin[nearx][iboot],xin[nearx+1][iboot]};
	    double tyin[3]={yin[nearx-1][iboot],yin[nearx][iboot],yin[nearx+1][iboot]};
	    
	    //find the spline
	    double a,b,c;
	    parabolic_spline(a,b,c,txin,tyin);
	    
	    //interpolate
	    yout[ixout].data[iboot]=a*sqr(xout[ixout][iboot])+b*xout[ixout][iboot]+c;

	    if(outpath!=NULL && ixout==0)
	      for(int iplot=0;iplot<nplot;iplot++)
		{
		  xplot[iplot]=xin[nearx-1][iboot]+(xin[nearx+1][iboot]-xin[nearx-1][iboot])/99*iplot;
		  if(iboot==0) cout<<xplot[iplot]<<endl;
		  yplot[iplot].data[iboot]=a*sqr(xplot[iplot])+b*xplot[iplot]+c;
		}
	  }
    }
  
  if(outpath!=NULL)
    {
      ofstream outplot(outpath);
      outplot<<"@type xydy"<<endl;
      outplot<<"@s0 symbol 1"<<endl;
      for(int ixin=0;ixin<nxin;ixin++)
	outplot<<xin[ixin].med()<<" "<<yin[ixin]<<endl;
      outplot<<"&\n@type xydxdy"<<endl;
      for(int ixout=0;ixout<nxout;ixout++)
	{
	  outplot<<xout[ixout].med()<<" "<<yout[ixout].med()<<" ";
	  outplot<<xout[ixout].err()<<" "<<yout[ixout].err()<<endl;
	}
      outplot<<"&"<<endl;
      write_polygon(outplot,xplot,yplot,1);
      for(int ixout=0;ixout<nxout;ixout++)
	{
	  outplot<<"&"<<endl;
	  for(int iboot=0;iboot<=nboot;iboot++)
	    outplot<<xout[ixout][iboot]<<" "<<yout[ixout][iboot]<<endl;
	}

      outplot.close();
    }
  
  return yout;
}

bvec interpolate_multi(double *xin,bvec yin,bvec xout,const char *outpath=NULL)
{
  bvec xin_vec(yin.nel,yin.nboot,yin.njack);
  for(int iel=0;iel<yin.nel;iel++)
    for(int iboot=0;iboot<=yin.nboot;iboot++)
      xin_vec[iel].data[iboot]=xin[iel];
  
  return interpolate_multi(xin_vec,yin,xout,outpath);
}

boot interpolate_single(bvec xin,bvec yin,boot xout,const char *outpath=NULL)
{
  bvec xout_vec(1,yin.nboot,yin.njack);
  xout_vec=xout;
  
  return interpolate_multi(xin,yin,xout_vec,outpath)[0];
}

boot interpolate_single(double *xin,bvec yin,boot xout,const char *outpath=NULL)
{
  bvec xout_vec(1,yin.nboot,yin.njack);
  xout_vec=xout;
  
  return interpolate_multi(xin,yin,xout_vec,outpath)[0];
}
