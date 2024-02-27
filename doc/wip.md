#### Paralleling transport whatever



``` c++

auto parallelTransport(const DerivedFromNode auto& conf,
                       const DerivedFromNode auto& v,
					   const Ori& ori,
                       const Dir& dir)
{
  return conf[dir]*shift(v,ori,dir);
}
```

* `v` and `conf` can be any data structure with `LocSite` index
* `conf` must have a Lorentz index


### Computing the plaquette

Compute the plaquette of whatever (for example, an `SU(3)` gauge configuration)


```c++

template <DerivedFromTransposableComp C=Color>     /// Support any Index
double plaquette(const DerivedFromNode auto& conf) /// Support any expression
{
  ComplexScalarField squares(0.0);
  
  for(Dir mu=0;mu<NDIM;mu++)
    for(Dir nu=mu+1;nu<NDIM;nu++)
      squares+=real(traceOver<C>((parallelTransport(conf,conf[nu],bw,mu)).close()*
				                 (parallelTransport(conf,conf[mu],bw,nu)).close()));
  
  const double plaq=
    squares.glbReduce()/lat->getGlbVol()/(2*sqr(C::size));
  
  return plaq;
}

```

* `conf` allocated on `gpu` -> kernels are issued
* `conf` allocated on `cpu` -> thread loop is issued
* Index order irrelevant
* Closing partial expression could be automatized


