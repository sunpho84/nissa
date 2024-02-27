### Tensor algebra

```c++
  DynamicTens<OfComps<LocLxSite,Dir>,double> tens(lat->getLocVol());`
```

allocates a tensor with a static-sized Lorentz index `Dir` and a dynamical-size specified index `LocLxSite`

One can add as many component as needed, and create new components

``` c++
  DECLARE_TRANSPOSABLE_COMP(Color,int,NCOL,color);
  DynamicTens<OfComps<LocLxSite,Dir,Color>,double> tens1(lat->getLocVol());
  DynamicTens<OfComps<Dir,LocLxSite,Color>,double> tens2(lat->getLocVol());
```

automatic understand tensor algebra

```c++
  DynamicTens<OfComps<LocLxSite,Color,Dir>,double> tens3=tens1*tens2;                            /// Direct product
  DynamicTens<OfComps<LocLxSite,ColorRow,ColorCln,DirRow,DirCln>,double> tens4=tens1*dag(tens2); /// Outer product
  DynamicTens<OfComps<LocLxSite>,double> tens5=dag(tens1)*tens2;                                 /// Scalar prodcut
```

Indices can be accessed via member methods or subscribing. Partial subscribing is supported, order of subscription is irrelevant. The constructions

```c++
  tens4.colorRow(0).locLxSite(52);
  tens4[ColorRow(0)][LocLxSite(52)];
  tens4(colorRow(0),locLxSite(52));
  tens4.locLxSite(52).colorRow(0);
```

are all equivalent, methods are automatically after the index is specified. Tensor support arbitrary underlying data type

```c++
  Struct FlavProperties
  {
     double charge;
  };

  DECLARE_UNTRANSPOSABLE_COMP(Flavor,int,NFLAV,flavor);
  DynamicTens<OfComps<FLAVOR>,FlavProperties> families;
  families.flavor(0).charge=0.3;

```

### Field

``` c++
  Field<OfComps<Dir,ColorRow,ColorCln,ComplId>,double> conf;
     
```

`conf` is a field with Lorentz index `Dir`, two color indices, a complex real/imag part, of basic type double, and implicit `LoclxSite` (spacetime)


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


