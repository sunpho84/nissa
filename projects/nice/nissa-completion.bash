# Bash completion for NICE; source this file to enable it
# Override NISSA_COMPLETION_PROGRAM when the query executable has another path

# Decode literal shell arguments up to point without expansion or evaluation
# Results use caller-local variables through Bash's dynamic scoping
_nissa_shell_prefix()
{
    local line=$1 c next quote= started=0 i valid=1
    _nissa_words=()
    _nissa_text= _nissa_quote= _nissa_drop=0
    for ((i=0; i<${#line} && valid; i++)); do
        c=${line:i:1}
        next=${line:i+1:1}
        if [[ $quote == "'" ]]; then
            if [[ $c == "'" ]]; then quote=; else _nissa_text+=$c; fi
        elif [[ $quote == '"' ]]; then
            case $c in
                '"') quote= ;;
                '$'|'`') valid=0 ;;
                '\')
                    if [[ $next == [\$\`\"\\] ]]; then
                        _nissa_text+=$next
                        ((i++))
                    elif [[ $next == $'\n' ]]; then
                        ((i++))
                    else
                        _nissa_text+=$c
                    fi ;;
                *) _nissa_text+=$c ;;
            esac
        else
            case $c in
                ' '|$'\t'|$'\n')
                    if ((started)); then _nissa_words+=("$_nissa_text"); fi
                    _nissa_text= _nissa_quote= _nissa_drop=0 started=0 ;;
                "'"|'"')
                    if ((${#_nissa_text} == _nissa_drop)) && [[ -z $_nissa_quote ]]; then
                        _nissa_quote=$c
                    else
                        valid=0
                    fi
                    quote=$c started=1 ;;
                '\')
                    if [[ -z $next ]]; then
                        valid=0
                    elif [[ $next != $'\n' ]]; then
                        _nissa_text+=$next
                        started=1
                    fi
                    ((i++)) ;;
                '$'|'`'|';'|'|'|'&'|'<'|'>'|'('|')') valid=0 ;;
                '#')
                    if ((started)); then _nissa_text+=$c; else valid=0; fi ;;
                *)
                    _nissa_text+=$c started=1
                    if [[ ${COMP_WORDBREAKS-} == *"$c"* ]]; then
                        _nissa_drop=${#_nissa_text}
                    fi ;;
            esac
        fi
    done
    if [[ -n $_nissa_quote && $quote != "$_nissa_quote" ]]; then valid=0; fi
    ((valid))
}

# Quote a replacement for the current Readline quote context
_nissa_quote_candidate()
{
    local value=$1
    if [[ $_nissa_quote == "'" ]]; then
        _nissa_candidate=${value//\'/\'\\\'\'}
    elif [[ $_nissa_quote == '"' ]]; then
        value=${value//\\/\\\\}
        value=${value//\"/\\\"}
        value=${value//\$/\\\$}
        _nissa_candidate=${value//\`/\\\`}
    else
        printf -v _nissa_candidate '%q' "$value"
    fi
}

_nissa_complete()
{
    # Bash reports character positions in its active locale; split before byte arithmetic
    local prefix=${COMP_LINE:0:COMP_POINT} after=${COMP_LINE:COMP_POINT}
    local LC_ALL=C
    local -a _nissa_words response
    local _nissa_text _nissa_quote _nissa_drop _nissa_candidate
    local argument context= options=1 file=0 i
    local program=${NISSA_COMPLETION_PROGRAM:-$1} request offset begin tail=
    COMPREPLY=()
    if ! _nissa_shell_prefix "$prefix"; then return; fi

    # Completed arguments determine option/file context; script pieces stay ordered
    for ((i=1; i<${#_nissa_words[@]}; i++)); do
        argument=${_nissa_words[i]}
        if ((file)); then
            file=0
        elif ((options)) && [[ $argument == -f ]]; then
            file=1
        elif ((options)) && [[ $argument == -- ]]; then
            options=0
        else
            context+="$argument"$'\n'
        fi
    done

    if ((file)); then
        # A quoted prefix keeps wildcard characters literal; glob results stay in array elements
        for argument in "$_nissa_text"*; do
            if [[ -e $argument || -L $argument ]]; then
                if [[ -d $argument ]]; then
                    argument+=/
                    compopt -o nospace 2>/dev/null || :
                fi
                COMPREPLY+=("${argument:_nissa_drop}")
            fi
        done
        compopt -o filenames 2>/dev/null || :
    elif ((options)) && [[ $_nissa_text == -* ]]; then
        for argument in -f -- --complete-at; do
            if [[ $argument == "$_nissa_text"* ]]; then
                _nissa_quote_candidate "$argument"
                COMPREPLY+=("$_nissa_candidate")
            fi
        done
    else
        # Preserve an existing identifier suffix when point is in the middle of a word
        while [[ $after == [a-zA-Z0-9_]* ]]; do
            tail+=${after:0:1}
            after=${after:1}
        done
        request=$context$_nissa_text$tail
        offset=$((${#context}+${#_nissa_text}))
        mapfile -d '' -t response < <(
            printf '%s' "$request" |
                "$program" --complete-at "$offset" --completion-format bash 2>/dev/null
        )
        if ((${#response[@]} >= 2)) && [[ ${response[0]} =~ ^[0-9]+$ && ${response[1]} =~ ^[0-9]+$ ]]; then
            begin=$((${response[0]}-${#context}))
            if ((begin >= _nissa_drop && begin <= ${#_nissa_text})); then
                for argument in "${response[@]:2}"; do
                    if [[ -z $tail || $argument == *"$tail" ]]; then
                        if [[ -n $tail ]]; then argument=${argument:0:${#argument}-${#tail}}; fi
                        _nissa_quote_candidate "${_nissa_text:_nissa_drop:begin-_nissa_drop}$argument"
                        COMPREPLY+=("$_nissa_candidate")
                    fi
                done
            fi
        fi
        compopt -o nospace -o noquote 2>/dev/null || :
    fi
}

complete -F _nissa_complete NICE
