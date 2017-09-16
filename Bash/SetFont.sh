#!/bin/bash

Font="Deja Vu Sans Mono"; Size=9; Style="Regular"

while getopts :f:s:l:h Option; do
    case "${Option}" in
        f)
            Font=${OPTARG}
        ;;
        s)
            Size=${OPTARG}
        ;;
        l)
            Style=${OPTARG}
        ;;
        h) 
            echo -en '\n'
            echo "Usage: $0 -f [FontName] | -s [FontSize] | -l [FontStyle]"
            echo -en '\n'
            
            exit 0
        ;;
        :)
            echo "Option -"$OPTARG" needs a value. See "$0 -h"."
            exit 0
        ;;
        \?)
            echo "Unknown option "$OPTARG". See "$0 -h"."
            exit 0
        ;;
    esac
done

printf '\33]50;%s\007' "xft:$Font:size=$Size:style=$Style"
