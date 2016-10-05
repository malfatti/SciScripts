#!/bin/bash

## Remove some special characters and change <> {} [] () to _
for i in *; do 
	new=`echo "$i" | \
		 tr -d '\"' | \
		 tr -d '\!' | \
		 tr -d '\@' | \
		 tr -d '\#' | \
		 tr -d '\$' | \
		 tr -d '\%' | \
		 tr -d '\&' | \
		 tr -d '\*' | \
		 tr -d '\+' | \
		 tr -d '\=' | \
		 tr -d '\,' | \
		 tr -d '\;' | \
		 tr -d '\:' | \
		 tr -d '\\' | \
		 tr -d '\/' | \
		 tr -d '\?' | \
		 tr -d '\€' | \
		 tr -d '\§' | \
		 tr -d '\¢' | \
		 tr -d '\£' | \
		 tr 'ç' c | \
		 tr 'Ç' C | \
		 tr '(' _ | \
		 tr ')' _ | \
		 tr '[' _ | \
		 tr ']' _ | \
		 tr '{' _ | \
		 tr '}' _ | \
		 tr '<' _ | \
		 tr '>' _`
	mv "$i" "$new"
done

## Remove single quotes and brackets
for i in "$(find . -type f -name "*'*")"; do 
	new=`echo "$i" | sed 's/\x27//g'`
	mv "$i" "$new"
done

## Lowercase file names
for i in *; do 
	mv "$i" "$(echo $i|tr A-Z a-z)"
done

## Capitalize every word of file names
for i in *; do 
	new=`echo "$i" | sed -e 's/^./\U&/g; s/ ./\U&/g'`
	mv "$i" "$new"
done

## Remove spaces on filenames
for i in *; do 
	new=`echo "$i" | sed 's/ //g'`
	mv "$i" "$new"
done
