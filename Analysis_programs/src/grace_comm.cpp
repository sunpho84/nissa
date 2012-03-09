grace grace::print_graph(double x,TYPE y)
{
  fout<<"@type xydy"<<endl;
  fout<<x<<" "<<y<<endl;

  return *this;
}

grace grace::print_graph(VTYPE y)
{
  int nel=y.nel;

  fout<<"@type xydy"<<endl;
  for(int i=0;i<nel;i++) fout<<i<<" "<<y[i]<<endl;

  return *this;
}

grace grace::print_graph(double *x,VTYPE y)
{
  int nel=y.nel;

  fout<<"@type xydy"<<endl;
  for(int i=0;i<nel;i++) fout<<x[i]<<" "<<y[i]<<endl;

  return *this;
}

void build_graph(double *x,VTYPE &y,double (*fun)(double,double*),double x0,double x1,int nx,VTYPE &par)
{
  double dx=(x1-x0)/nx;
#ifdef JACK
  VTYPE xj(nx+1,par.njack);
  y=VTYPE(nx+1,par.njack);
#else
  VTYPE xj(nx+1,par.nboot,par.njack);
  y=VTYPE(nx+1,par.nboot,par.njack);
#endif
  
  for(int i=0;i<=nx;i++) xj.data[i]=x[i]=x0+dx*i;
  
  y=par_single_fun(fun,xj,par);
} 

///////////////////////////////////////////////////////////////////////////////

void write_polygon(ofstream &fout,double *x,VTYPE &in,int nset=0)
{
  int nel=in.nel,N=in.N;
  double *err=new double[nel];

  fout<<"@type xy"<<endl;
  fout<<"@s"<<nset<<" fill type 1"<<endl;
  
  //lower line
  for(int i=0;i<nel;i++)
    {
      err[i]=in[i].err();
      fout<<x[i]<<" "<<in.data[i].data[N]-err[i]<<endl;
    }
  
  //upper line
  for(int i=nel-1;i>=0;i--)fout<<x[i]<<" "<<in.data[i].data[N]+err[i]<<endl;
  
  delete[] err;
}

grace grace::polygon(double *x,VTYPE &in)
{
  write_polygon(fout,x,in,nset);
  return *this;
}

grace grace::polygon(double (*fun)(double,double*),double x0,double x1,int nx,VTYPE &par)
{
  double x[nx];
  VTYPE y;
  
  build_graph(x,y,fun,x0,x1,nx,par);
  polygon(x,y);
  
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

grace grace::contour(double *x,VTYPE &in)
{
  int nel=in.nel,N=in.N;
  double *err=new double[nel];

  fout<<"@type xy"<<endl;

  //lower line
  for(int i=0;i<nel;i++)
    {
      err[i]=in[i].err();
      fout<<x[i]<<" "<<in.data[i].data[N]-err[i]<<endl;
    }

  //upper line  
  for(int i=nel-1;i>=0;i--) fout<<x[i]<<" "<<in.data[i].data[N]+err[i]<<endl;
  
  //junction
  fout<<x[0]<<" "<<in.data[0].data[N]-err[0]<<endl;

  //central line
  for(int i=0;i<nel;i++) fout<<x[i]<<" "<<in.data[i].data[N]<<endl;

  
  delete[] err;
  return *this;
}

grace grace::contour(double (*fun)(double,double*),double x0,double x1,int nx,VTYPE &par)
{
  double x[nx];
  VTYPE y;
  
  build_graph(x,y,fun,x0,x1,nx,par);
  contour(x,y);
  
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////

grace grace::ave_line(double *x,VTYPE &in)
{
  int nel=in.nel;
  double y[nel];
  
  for(int iel=0;iel<nel;iel++) y[iel]=in.data[iel].med();
  
  return ave_line(x,y,nel);
}

