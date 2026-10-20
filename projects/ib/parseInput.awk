BEGIN{
    print "digraph test {"
}

{
    if($3=="LINCOMB")
    {
	offs=1+$4*2
	source="LINCOMB"NR
	for(i=5;i<4+$4*2;i+=2)
	    print $i"->LINCOMB"NR" [label=\""$(i+1)"\"]"
    }
    else
    {
	label=$2
	offs=0
	source=$3
    }
}

{
    t=$(offs+4)
    k=$(offs+5)
    m=$(offs+6)
    r=$(offs+7)
    
    maybeT=(t==-1)?"":",t="t
}

/* Identify props */
$2=="-"{
    label="PROP(k="k",m="m",r="r""maybeT")"
}

/* Identify gamma insertion */
$2=="G"{
    label="G("k""maybeT")"
}

/* Identify smearing */
$2=="SM"{
    kl=""
    kx=$(offs+7)
    if(kx)
	kl=kl",px="kx
    ky=$(offs+8)
    if(ky)
	kl=kl",py="ky
    kz=$(offs+9)
    if(kz)
	kl=kl",pz="kz
    label="SM(k="k",r="m""maybeT""kl")"
}

$1!="/*" && NF>1 && n && b==2{
    print source,"->",$1,"[label=\""label"\"]"
    n--
}

/* Capture the number of props */
$1=="NProps"{
    n=$2
}

/* Capture the beginning of props: needs to get two "Name" as the first is the source list */
$1=="Name"{
    b++
}
nm>0 && NF==3{
    print $2,"->",$1""nm
    print $3,"->",$1""nm
    nm--
}

$1=="NMes2PtsContr"{
    nm=$2
}

END{
    print "}"
}
