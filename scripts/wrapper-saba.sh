#!/bin/bash

PIPELINE_HOME=/home/dresdnerg/pipeline
R=/opt/bin/R

if [ -z "$SGE_TASK_ID" ];
then
        SGE_TASK_ID="0"
    fi

    i=1
    for f in `ls /home/dresdnerg/pipeline/cbs/in/*.cleaned | sort`;
    do
            if [ $i -eq $SGE_TASK_ID ];
                    then
                                #/usr/bin/python $PIPELINE_HOME/scripts/print_mapped.py $PIPELINE_HOME/markerfiles/cleaned_markers.txt $PIPELINE_HOME/normalized/$f /home/dresdnerg/pipeline/cbs/in/
                                        $R CMD BATCH --no-save "--args $f $PIPELINE_HOME/cbs/out/" /home/dresdnerg/pipeline/scripts/run_cbs.r $f.cbs.log
                                            fi
                                                let i=i+1;
                                            done
            fi
    done
fi
