#pragma once

typedef char gr_string[20];

class grace
{
public:
  int nset;
  ofstream fout;
  ostream& operator<<(jvec out);
  ostream& operator<<(bvec out);
  
  grace(const char *path,...);
  grace(const grace&);
  grace();
  
  void close(){fout.close();}
  virtual ~grace(){close();}
  
  grace plot_size(int,int);
  grace plot_title(const char*);
  
  grace xaxis_label(const char*);
  grace yaxis_label(const char*);
  grace axis_label(const char*,const char*);

  grace set_line_size(int);  
  grace set_line_type(const char*);
  grace set_color(const char*);
  grace set_legend(const char*);
  grace set_symbol(const char*);
  grace new_set(...);
  grace set(int,...);
  
  grace ave_line(double *,double *,int);
  
  grace print_graph(double *,jvec);
  grace print_graph(double,jack);
  grace print_graph(jvec);

  grace contour(double *,jvec &);
  grace polygon(double *,jvec &);
  grace ave_line(double *,jvec &);
  grace polygon(double (*fun)(double,double*),double,double,int,jvec &);
  grace contour(double (*fun)(double,double*),double,double,int,jvec &);

  grace print_graph(double *,bvec);
  grace print_graph(double,boot);
  grace print_graph(bvec);

  grace contour(double *,bvec &);
  grace polygon(double *,bvec &);
  grace ave_line(double *,bvec &);
  grace contour(double (*fun)(double,double*),double,double,int,bvec &);
  grace polygon(double (*fun)(double,double*),double,double,int,bvec &);
};

int switch_gr_string(const gr_string *list,int nel,const char *in)
{
  int ind=0;
  while(strcmp(list[ind],in)&&ind<nel) ind++;
  if(ind==nel) ind=-1;
  
  return ind;
}

int switch_color(const char* col)
{
  const int ngr_col=16;
  const gr_string gr_col[ngr_col]={"white","black","red","green","blue","yellow","brown","grey","violet","cyan","magenta","orange","indigo","maroon","turquoise","green4"};

  return switch_gr_string(gr_col,ngr_col,col);
}

int switch_symbol(const char* symb)
{
  const int ngr_symbol=11;
  const gr_string gr_symb[ngr_symbol]={"none","circle","square","diamond","triangle up","triangle left","triangle down","plus","x","star","char"};
  
  return switch_gr_string(gr_symb,ngr_symbol,symb);
}

int switch_line_type(const char* type)
{
  const int ngr_line_type=3;
  const gr_string gr_line_type[ngr_line_type]={"straight","pointed","dashed"};

  return switch_gr_string(gr_line_type,ngr_line_type,type);
}

grace::grace(const grace &in) : nset(in.nset)
{fout.tie(in.fout.tie());}

grace::grace()
{
  cerr<<"Error, cannot create a grace file untied to anything!"<<endl;
  exit(1);
}

grace::grace(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  nset=0;
  fout.open(buffer);
  if(!fout.good())
    {
      cerr<<"Errorr in opening grace file: "<<buffer<<endl;
      exit(1);
    }
}

ostream& grace::operator<<(jvec out)
{return fout<<out;}

ostream& grace::operator<<(bvec out)
{return fout<<out;}

grace grace::new_set(...)
{
  fout<<"&"<<endl;
  nset++;
  return *this;
}

grace grace::ave_line(double *x,double *y,int n)
{
  fout<<"@type xy"<<endl;
  for(int i=0;i<n;i++) fout<<x[i]<<" "<<y[i]<<endl;
  return *this;
}

#define JACK
#define VTYPE jvec
#define TYPE jack
#define N njack
#include "grace_comm.cpp"
#undef N
#undef TYPE
#undef VTYPE
#undef JACK

#define BOOT
#define VTYPE bvec
#define TYPE boot
#define N nboot
#include "grace_comm.cpp"
#undef N
#undef TYPE
#undef VTYPE
#undef BOOT

grace grace::plot_size(int x,int y)
{
  fout<<"@page size "<<x<<", "<<y<<endl;
  return *this;
}

grace grace::plot_title(const char *title)
{
  fout<<"@title \""<<title<<"\""<<endl;
  return *this;
}

grace grace::xaxis_label(const char *title)
{
  fout<<"@xaxis  label \""<<title<<"\""<<endl;
  return *this;
}

grace grace::yaxis_label(const char *title)
{
  fout<<"@yaxis  label \""<<title<<"\""<<endl;
  return *this;
}

grace grace::axis_label(const char *titlex,const char *titley)
{
  fout<<"@xaxis  label \""<<titlex<<"\""<<endl;
  fout<<"@yaxis  label \""<<titley<<"\""<<endl;
  return *this;
}

grace grace::set_legend(const char *leg)
{
  fout<<"@s"<<nset<<" legend \""<<leg<<"\""<<endl;
  return *this;
}

grace grace::set_color(const char *col)
{
  int icol=switch_color(col);

  if(icol!=-1)
    {
      fout<<"@s"<<nset<<" symbol color "<<icol<<endl;
      fout<<"@s"<<nset<<" symbol fill color "<<icol<<endl;
      
      fout<<"@s"<<nset<<" line color "<<icol<<endl;
      fout<<"@s"<<nset<<" fill color "<<icol<<endl;
      
      fout<<"@s"<<nset<<" avalue color "<<icol<<endl;
      fout<<"@s"<<nset<<" errorbar color "<<icol<<endl;
    }
  
  return *this;
}

grace grace::set_symbol(const char *symb)
{
  int type=switch_symbol(symb);
  if(type!=-1) fout<<"@s"<<nset<<" symbol "<<type<<endl;
  
  if(strcmp(symb,"filled")==0) fout<<"@s"<<nset<<" symbol fill pattern 1"<<endl;
  if(strcmp(symb,"unfilled")==0) fout<<"@s"<<nset<<" symbol fill pattern 0"<<endl;
  
  return *this;
}

grace grace::set_line_size(int line_size)
{
  fout<<"@s"<<nset<<" symbol linewidth "<<line_size<<endl;  
  fout<<"@s"<<nset<<" line linewidth "<<line_size<<endl;  
  fout<<"@s"<<nset<<" errorbar linewidth "<<line_size<<endl;  
  fout<<"@s"<<nset<<" errorbar riser linewidth "<<line_size<<endl;  
  return *this;
}

grace grace::set_line_type(const char *line_type)
{
  if(!strcmp(line_type,"none")) fout<<"@s"<<nset<<" line type 0"<<endl;
  else
    {
      int type=switch_line_type(line_type);
      if(type!=-1)
	{
	  fout<<"@s"<<nset<<" line type 1"<<endl;
	  fout<<"@s"<<nset<<" line type "<<type;
	}
    }
  
  return *this;
}

grace grace::set(int n,...)
{
  char *what;

  va_list vl;
  va_start(vl,n);
  for(int i=0;i<n;i++)
    {
      what=va_arg(vl,char*);
      set_color(what);
      set_line_type(what);
      set_symbol(what);
    }
  va_end(vl);
  
  return *this;
}
