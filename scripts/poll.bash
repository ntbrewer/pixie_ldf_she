# /bin/bash
if [[ -z $1 ]];then
    echo "Must provide at least one argument"
    exit 1;
fi
 
if [[ ${1:0:1} != - ]];then
   echo "Usage error."
   echo $1 "needs leading - for options"
   exit 1;
fi

if [[ ${1:(-1)} == t ]];then
    if [[ -z $2 || $2 -le 0  ]];then
        echo "Usage error."
        echo "Use number >0 after option -t"
        exit 1;
    fi
fi

echo "ok."
var=1




until ./poll $1 $2; do
     SUB="Poll Crashed!!!!"
     EmailAddr="brewer.nathant@gmail.com"
     EmailText="/tmp/pollemail.txt"
     echo "Poll crashed and is trying to restart." > $EmailText
     echo "Check it! And find out why!" >> $EmailText
     /bin/mail -s "$SUB" "$EmailAddr" < $EmailText
    sleep 1
    ((var++))
    echo $var
    if [[ $var == 3 ]];then
         SUB="Poll Crashed!!!!"
         EmailAddr="brewer.nathant@gmail.com, 6157123925@messaging.sprintpcs.com"
         EmailText="/tmp/pollemail.txt"
         echo "Poll crashed and cannot restart." > $EmailText
         echo "Check it! And find out why!" >> $EmailText
         /bin/mail -s "$SUB" "$EmailAddr" < $EmailText
         exit 1;
    fi
done
