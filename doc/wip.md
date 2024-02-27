### Computing the plaquette

Compute the plaquette of whatever (for example, an `SU(3)` gauge configuration)

```c++

template <DerivedFromTransposableComp C=Color>     /// Support any Index
double plaquette(const DerivedFromNode auto& conf) /// Support any expression
{
  ComplexScalarField squares(0.0);
  
  for(Dir mu=0;mu<NDIM;mu++)
    for(Dir nu=mu+1;nu<NDIM;nu++)
      squares+=real(traceOver<C>((conf[mu]*shift(conf[nu],bw,mu)).close()*
                                 (conf[nu]*shift(conf[mu],bw,nu)).close()));
  
  const double plaq=
    squares.glbReduce()/lat->getGlbVol()/(2*sqr(C::size));
  
  return plaq;
}

```

* Closing partial expression could be automatized
* Global reduction done in a trivial way
