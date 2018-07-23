files=`ls`
for fname in $files ; do mv "$fname" "$(echo "$fname" | sed s/0.//g)" ; done
