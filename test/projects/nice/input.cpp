{
  struct Smear
  {
    kappa=4.0;
    
    n;
    
    conf;
  };
  
  fun smear(n)
  {
    return Smear(.n=n,.conf="ciao");
  }
  
  struct Source
  {
    spinDiluted;
    
    colorDiluted;
    
    noiseType;
    
    id;
  };
  
  fun source(id)
  {
    return Source(.spinDiluted=1,.colorDiluted=0,.noiseType=0,.id=id);
  }
  
  struct QProp
  {
    kappa;
    
    mu;
    
    residue;
    
    csw;
    
    thetaX=0.0;
    
    thetaY=0.0;
    
    thetaZ=0.0;
  };
  
  glbCsw=1.69;
  
  glbKappa=0.1214;
  
  fun qProp(mu,residue=1e-16,thetaX=0.0,thetaY=0.0,thetaZ=0.0)
  {
    return QProp(.kappa=glbKappa,.mu=mu,.residue=residue,.csw=glbCsw,.thetaX=thetaX,.thetaY=thetaY,.thetaZ=thetaZ);
  }
  
  struct SelectT
  {
    t;
  };
  
  struct Prod
  {
    a;
    
    b;
  };
  
  fun prod(a,b)
  {
    return Prod(.a=a,.b=b);
  }
  
  SM10=smear(10);
  
  eta=source(0);
  
  L=qProp(0.0123);
  
  print(L,"\n");
  
  line=prod(prod(L,SM10),eta);
  
  print(line,"\n");
  
  er=1;
  er=2;
  1+2-3/2*2%5<3>1<=4>=6==3!=8;
  !2.0;
  er++;
  er--;

  fun cicc(ar,&yt,vf=1,arrrrg=1)
  {
    yt[2]=12;
  }
  
  l=lambda(){};
  v=lambda(...){};
  
  s=seq(10);
  
  a=s[9];
  print(a,"\n");
  
  a*=3;
  
  p=lambda(x)
    {
      print(x+1,"\n");
    };
  p(a);
  forEachEl(s,p);
  cicc(1,&s);
  print(s[2],"\n");
  
  struct Ss
  {
    d=1;
    
    e;
  };
  
  ss=Ss(.e=2);
  
  print("ss.d: ",ss.d,"\n");
  
  ss.d=2;
  
  print("ss.d: ",ss.d,"\n");
  print("ss.e: ",ss.e,"\n");
  
  tt=Ss(.e=9);
  
  print("tt.d: ",tt.d,"\n");
  print("tt.e: ",tt.e,"\n");
  print("ss.e: ",ss.e,"\n");
}
