#!/bin/csh

rm -f CPUTIME.txt

set ERRS = (1e-2 1e-3 1e-4 1e-5)
set DIRS = (mmex4 mmex6 mmex8 mmex10 fdex4 fdex6 fdex8 fdex10)
set STATS = (Mean Min Max SD)

foreach STAT ($STATS)

    echo "$STAT" >> CPUTIME.txt

    # Column headings
    printf "%-8s" "Error" >> CPUTIME.txt
    foreach dir ($DIRS)
        printf " %12s" "$dir" >> CPUTIME.txt
    end
    printf "\n" >> CPUTIME.txt

    foreach ERR ($ERRS)

        printf "%-8s" "$ERR" >> CPUTIME.txt

        foreach dir ($DIRS)

            set value = `grep Sol ${dir}/run*${ERR}* | awk -v stat="$STAT" '{x=$8; sum+=x; sumsq+=x*x; if(n==0 || x<min) min=x; if(n==0 || x>max) max=x; n++} END{if(n==0){print "NaN"; exit} mean=sum/n; var=sumsq/n-mean*mean; if(var<0 && var>-1.0e-14) var=0; sd=sqrt(var); if(stat=="Mean") printf "%.7g",mean; else if(stat=="Min") printf "%.7g",min; else if(stat=="Max") printf "%.7g",max; else if(stat=="SD") printf "%.7g",sd}'`

            printf " %12s" "$value" >> CPUTIME.txt

        end

        printf "\n" >> CPUTIME.txt
    end

    printf "\n" >> CPUTIME.txt
end
