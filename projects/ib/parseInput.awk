BEGIN{
    print "digraph test {"
}

{
    label=$2
}

/* Identify props */
$2=="-"{
    label="PROP(k="$5",m="$6",r="$7")"
}

/* Identify gamma insertion */
$2=="G"{
label="G"$5}

/* Identify smearing */
$2=="SM"{
    label="SM(k="$5",m="$6")"
}

$1!="/*" && NF>1 && n && b==2{
    print $3,"->",$1,"[label=\""label"\"]"
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

END{
    print "}"
}
