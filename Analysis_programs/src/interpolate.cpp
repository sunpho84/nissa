#pragma once

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

boot interpolate_single(bvec vec,double *x,boot xf)
{
  int nx=vec.nel;
  int nboot=vec.nboot;
  int njack=vec.njack;
  boot out(nboot,njack);

  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      //check if the needed x is not below the lower one or above the heavier one
      double dmin=(x[1]-x[0])/2;
      double dmax=(x[nx-1]-x[nx-2])/2;
      bool linbelow=(xf[iboot]-x[0]<dmin);
      bool linabove=(x[nx-1]-xf[iboot]<dmax);
      if(nx==2) linabove=1;
      
      //if linear
      if(linabove||linbelow)
	{
	  //choose the pair
	  int sup;
	  if(linabove) sup=nx-1;
	  else sup=1;
	  
	  //find the m and q
	  double m=(vec[sup][iboot]-vec[sup-1][iboot])/(x[sup]-x[sup-1]);
	  double q=vec[sup-1][iboot]-m*x[sup-1];
	  
	  //interpolate
	  out.data[iboot]=m*xf[iboot]+q;
	  
	  if(iboot==nboot && 0)
	    {
	      cout<<" "<<sup-1<<" "<<x[sup-1]<<" "<<vec[sup-1][iboot]<<" "<<m*x[sup-1]+q<<endl;
	      cout<<" "<<sup<<" "<<x[sup]<<" "<<vec[sup][iboot]<<" "<<m*x[sup]+q<<endl;
	      cout<<"  int: "<<xf[iboot]<<" "<<out.data[iboot]<<endl;
	    }	  
	}
      else
	{
	  //if parabolic interpolation, find the nearest point
	  int nearx=0;
	  for(int ix=0;ix<nx;ix++) if(fabs(x[ix]-xf[iboot])<fabs(x[nearx]-xf[iboot])) nearx=ix;
	  
	  //copy the x and y for current bootstrap
	  double xin[3]={x[nearx-1],x[nearx],x[nearx+1]};
	  double yin[3]={vec[nearx-1][iboot],vec[nearx][iboot],vec[nearx+1][iboot]};
	  
	  //find the spline
	  double a,b,c;
	  parabolic_spline(a,b,c,xin,yin);
	  
	  //interpolate
	  out.data[iboot]=a*sqr(xf[iboot])+b*xf[iboot]+c;

	  if(iboot==nboot && 0)
	    {
	      cout<<" "<<nearx-1<<" "<<x[nearx-1]<<" "<<vec[nearx-1][iboot]<<endl;
	      cout<<" "<<nearx<<" "<<x[nearx]<<" "<<vec[nearx][iboot]<<endl;
	      cout<<" "<<nearx+1<<" "<<x[nearx+1]<<" "<<vec[nearx+1][iboot]<<endl;
	      cout<<"  int: "<<xf[iboot]<<" "<<out.data[iboot]<<endl;
	    }	  
	}
    }
  
  return out;
}
